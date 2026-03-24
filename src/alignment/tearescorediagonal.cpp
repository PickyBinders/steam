#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "DistanceCalculator.h"
#include "QueryMatcher.h"
#include "FastSort.h"

#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

// Log-linear E-value model following Edgar & Sahakyan (2025).
// E(s) = (H/Q) * 10^(m*s + c), where s = raw * min(qcov, tcov).
// See evalue_calibration/17_sweep_features.py.
static const double LOGLINEAR_M_LN = -0.017364 * 2.302585093;
static const double LOGLINEAR_C_LN = -0.7636 * 2.302585093;

static double computeEvalue(double score, double hitsPerQuery) {
    return hitsPerQuery * exp(LOGLINEAR_M_LN * score + LOGLINEAR_C_LN);
}

template<typename T>
static DistanceCalculator::LocalAlignment ungappedAlignment(const T *seqTea1,
                                                            const T *seqAA1,
                                                            const T *seqTea2,
                                                            const T *seqAA2,
                                                            const unsigned int length,
                                                            short **subTeaMat,
                                                            short **subAAMat) {
    int maxScore = 0;
    int maxEndPos = 0;
    int maxStartPos = 0;
    int minPos = -1;
    int score = 0;
    for (unsigned int pos = 0; pos < length; pos++) {
        int currTea = subTeaMat[static_cast<int>(seqTea1[pos])][static_cast<int>(seqTea2[pos])];
        int currAA = subAAMat[static_cast<int>(seqAA1[pos])][static_cast<int>(seqAA2[pos])];
        score = currTea + currAA + score;
        const bool isMinScore = (score <= 0);
        score = (isMinScore) ? 0 : score;
        minPos = (isMinScore) ? pos : minPos;
        const bool isNewMaxScore = (score > maxScore);
        maxEndPos = (isNewMaxScore) ? pos : maxEndPos;
        maxStartPos = (isNewMaxScore) ? minPos + 1 : maxStartPos;
        maxScore = (isNewMaxScore) ? score : maxScore;
    }
    return DistanceCalculator::LocalAlignment(maxStartPos, maxEndPos, maxScore);
}

static Matcher::result_t ungappedAlignTea(Sequence &qSeqAA, Sequence &qSeqTea,
                                           Sequence &tSeqAA, Sequence &tSeqTea,
                                           int diagonal, SubstitutionMatrix &subMatAA,
                                           SubstitutionMatrix &subMatTea,
                                           double hitsPerQuery,
                                           std::string &backtrace, const Parameters &par) {
    DistanceCalculator::LocalAlignment res;
    float seqId = 0.0;
    backtrace.clear();
    unsigned int minDistToDiagonal = abs(diagonal);
    res.distToDiagonal = minDistToDiagonal;
    res.diagonal = diagonal;

    if (diagonal >= 0 && minDistToDiagonal < static_cast<unsigned int>(qSeqAA.L)) {
        unsigned int minSeqLen = std::min(static_cast<unsigned int>(tSeqAA.L),
                                          static_cast<unsigned int>(qSeqAA.L) - minDistToDiagonal);
        res.diagonalLen = minSeqLen;
        DistanceCalculator::LocalAlignment tmp = ungappedAlignment(
            qSeqTea.numSequence + minDistToDiagonal,
            qSeqAA.numSequence + minDistToDiagonal,
            tSeqTea.numSequence,
            tSeqAA.numSequence,
            minSeqLen, subMatTea.subMatrix, subMatAA.subMatrix);
        res.score = tmp.score;
        res.startPos = tmp.startPos;
        res.endPos = tmp.endPos;
    } else if (diagonal < 0 && minDistToDiagonal < static_cast<unsigned int>(tSeqAA.L)) {
        unsigned int minSeqLen = std::min(static_cast<unsigned int>(tSeqAA.L) - minDistToDiagonal,
                                          static_cast<unsigned int>(qSeqAA.L));
        res.diagonalLen = minSeqLen;
        DistanceCalculator::LocalAlignment tmp = ungappedAlignment(
            qSeqTea.numSequence, qSeqAA.numSequence,
            tSeqTea.numSequence + minDistToDiagonal,
            tSeqAA.numSequence + minDistToDiagonal,
            minSeqLen, subMatTea.subMatrix, subMatAA.subMatrix);
        res.score = tmp.score;
        res.startPos = tmp.startPos;
        res.endPos = tmp.endPos;
    }

    unsigned int distanceToDiagonal = res.distToDiagonal;

    int qStartPos, qEndPos, dbStartPos, dbEndPos;
    if (diagonal >= 0) {
        qStartPos = res.startPos + distanceToDiagonal;
        qEndPos = res.endPos + distanceToDiagonal;
        dbStartPos = res.startPos;
        dbEndPos = res.endPos;
    } else {
        qStartPos = res.startPos;
        qEndPos = res.endPos;
        dbStartPos = res.startPos + distanceToDiagonal;
        dbEndPos = res.endPos + distanceToDiagonal;
    }

    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, dbStartPos, dbEndPos);
    if (par.addBacktrace) {
        backtrace.append(alnLength, 'M');
    }

    float queryCov = (std::min((unsigned int)qSeqAA.L, (unsigned int)qEndPos) - (unsigned int)qStartPos + 1) / (float)qSeqAA.L;
    float targetCov = (std::min((unsigned int)tSeqAA.L, (unsigned int)dbEndPos) - (unsigned int)dbStartPos + 1) / (float)tSeqAA.L;

    // Coverage-weighted score for E-value (sqrt coverage)
    double covScore = res.score * sqrt(std::min(queryCov, targetCov));
    double evalue = computeEvalue(covScore, hitsPerQuery);

    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov));
    if (hasLowerCoverage) {
        return Matcher::result_t(UINT_MAX, res.score, queryCov, targetCov, seqId, evalue, alnLength,
                                 qStartPos, qEndPos, qSeqAA.L, dbStartPos, dbEndPos, tSeqAA.L, backtrace);
    }
    // No E-value filtering in ungapped rescoring — let gapped alignment decide.
    {
        int idCnt = 0;
        for (int i = qStartPos; i <= qEndPos; i++) {
            char qLetter = qSeqAA.numSequence[i];
            char tLetter = tSeqAA.numSequence[dbStartPos + (i - qStartPos)];
            idCnt += (qLetter == tLetter) ? 1 : 0;
        }
        seqId = Util::computeSeqId(par.seqIdMode, idCnt, qSeqAA.L, tSeqAA.L, alnLength);
    }
    return Matcher::result_t(tSeqAA.getDbKey(), res.score, queryCov, targetCov, seqId, evalue, alnLength,
                             qStartPos, qEndPos, qSeqAA.L, dbStartPos, dbEndPos, tSeqAA.L, backtrace);
}

int tearescorediagonal(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    bool sameDB = (par.db1.compare(par.db2) == 0);

    // Load query TEA database (main DB)
    IndexReader qTeaDbr(par.db1, par.threads, IndexReader::SEQUENCES,
                        touch ? IndexReader::PRELOAD_INDEX : 0);

    // Load query AA database (_aa suffix)
    std::string qAaDbName = par.db1 + "_aa";
    DBReader<unsigned int> qAaDbr(qAaDbName.c_str(), (qAaDbName + ".index").c_str(), par.threads,
                                   DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qAaDbr.open(DBReader<unsigned int>::NOSORT);

    // Load target databases
    IndexReader *tTeaDbr = NULL;
    DBReader<unsigned int> *tAaDbr = NULL;
    if (sameDB) {
        tTeaDbr = &qTeaDbr;
        tAaDbr = &qAaDbr;
    } else {
        tTeaDbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES,
                                   touch ? IndexReader::PRELOAD_INDEX : 0);
        std::string tAaDbName = par.db2 + "_aa";
        tAaDbr = new DBReader<unsigned int>(tAaDbName.c_str(), (tAaDbName + ".index").c_str(), par.threads,
                                             DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        tAaDbr->open(DBReader<unsigned int>::NOSORT);
    }

    // Load prefilter results
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                         DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads),
                 par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    // TEA substitution matrix (from --matcha)
    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "TEA substitution matrix (--matcha) is required for tearescorediagonal\n";
        EXIT(EXIT_FAILURE);
    }
    SubstitutionMatrix subMatTea(par.teaMatrixFile.c_str(), 1.0, par.scoreBias);

    // AA substitution matrix (from --sub-mat, weighted by --aa-weight)
    float aaFactor = par.teaWeight;
    SubstitutionMatrix subMatAA(par.scoringMatrixFile.values.aminoacid().c_str(), aaFactor, par.scoreBias);

    Debug::Progress progress(resultReader.getSize());

    double hitsPerQuery = static_cast<double>(par.maxResListLen);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        Sequence qSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *)&subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeqTea(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *)&subMatTea, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *)&subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeqTea(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *)&subMatTea, 0, false, par.compBiasCorrection);

        std::string backtrace;
        char buffer[1024 + 32768];
        std::string resultBuffer;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);

            if (*data != '\0') {
                unsigned int queryId = qTeaDbr.sequenceReader->getId(queryKey);

                char *querySeqTea = qTeaDbr.sequenceReader->getData(queryId, thread_idx);
                char *querySeqAA = qAaDbr.getData(qAaDbr.getId(queryKey), thread_idx);
                unsigned int querySeqLen = qTeaDbr.sequenceReader->getSeqLen(queryId);

                qSeqTea.mapSequence(id, queryKey, querySeqTea, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    hit_t prefHit = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = prefHit.seqId;
                    unsigned int targetId = tTeaDbr->sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB));

                    char *targetSeqTea = tTeaDbr->sequenceReader->getData(targetId, thread_idx);
                    char *targetSeqAA = tAaDbr->getData(tAaDbr->getId(dbKey), thread_idx);
                    const int targetSeqLen = static_cast<int>(tTeaDbr->sequenceReader->getSeqLen(targetId));

                    tSeqTea.mapSequence(targetId, dbKey, targetSeqTea, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);

                    if (Util::canBeCovered(par.covThr, par.covMode, static_cast<int>(querySeqLen), targetSeqLen) == false) {
                        rejected++;
                        continue;
                    }

                    Matcher::result_t res = ungappedAlignTea(qSeqAA, qSeqTea,
                                                              tSeqAA, tSeqTea,
                                                              static_cast<short>(prefHit.diagonal),
                                                              subMatAA, subMatTea, hitsPerQuery,
                                                              backtrace, par);

                    if (res.dbKey == UINT_MAX) {
                        rejected++;
                        continue;
                    }

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        alignmentResult.emplace_back(res);
                        passedNum++;
                        rejected = 0;
                    } else {
                        rejected++;
                    }
                }
            }

            if (alignmentResult.size() > 1) {
                SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
            }
            for (size_t result = 0; result < alignmentResult.size(); result++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[result], par.addBacktrace);
                resultBuffer.append(buffer, len);
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }
    }

    dbw.close();
    resultReader.close();
    qAaDbr.close();

    if (!sameDB) {
        delete tTeaDbr;
        tAaDbr->close();
        delete tAaDbr;
    }

    return EXIT_SUCCESS;
}
