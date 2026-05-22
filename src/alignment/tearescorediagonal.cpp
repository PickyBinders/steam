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

// Reuse mmseqs's parsePrecisionLib via extern (no `static` in upstream).
extern float parsePrecisionLib(const std::string &scoreFile, double targetSeqid, double targetCov, double targetPrecision);

// Pull in the precision-library data. Wrap in an anonymous namespace so the
// non-static `_lib_len` does not collide with mmseqs's own definition.
namespace {
#include "CovSeqidQscPercMinDiag.lib.h"
#include "CovSeqidQscPercMinDiagTargetCov.lib.h"
}

#ifdef OPENMP
#include <omp.h>
#endif

// Log-linear E-value model following Edgar & Sahakyan (2025).
// E(s) = P(FP) * (H/Q) * 10^(m*s + c), where s = raw alignment score.
// Parameters configurable via --loglinear-m, --loglinear-c, --p-fp.
static double computeEvalue(double rawScore, double hitsPerQuery, double m_ln, double c_ln, double pfp) {
    return pfp * hitsPerQuery * exp(m_ln * rawScore + c_ln);
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
                                           double hitsPerQuery, double m_ln, double c_ln, double pfp,
                                           std::string &backtrace, const Parameters &par,
                                           float &outScorePerCol) {
    DistanceCalculator::LocalAlignment res;
    float seqId = 0.0;
    backtrace.clear();
    unsigned int minDistToDiagonal = abs(diagonal);
    res.distToDiagonal = minDistToDiagonal;
    res.diagonal = diagonal;
    outScorePerCol = 0.0f;

    // HAMMING mode (used by linclust's pre-cluster step): no local SW,
    // just count AA identities along the entire diagonal overlap.
    if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
        unsigned int qOff = (diagonal >= 0) ? minDistToDiagonal : 0;
        unsigned int tOff = (diagonal >= 0) ? 0 : minDistToDiagonal;
        if (qOff >= static_cast<unsigned int>(qSeqAA.L) || tOff >= static_cast<unsigned int>(tSeqAA.L)) {
            return Matcher::result_t(UINT_MAX, 0, 0, 0, 0, 0, 0, 0, 0, qSeqAA.L, 0, 0, tSeqAA.L, backtrace);
        }
        unsigned int diagLen = std::min(static_cast<unsigned int>(qSeqAA.L) - qOff,
                                        static_cast<unsigned int>(tSeqAA.L) - tOff);
        int idCnt = 0;
        for (unsigned int i = 0; i < diagLen; i++) {
            idCnt += (qSeqAA.numSequence[qOff + i] == tSeqAA.numSequence[tOff + i]) ? 1 : 0;
        }
        seqId = Util::computeSeqId(par.seqIdMode, idCnt, qSeqAA.L, tSeqAA.L, diagLen);
        float queryCov = (float)diagLen / (float)qSeqAA.L;
        float targetCov = (float)diagLen / (float)tSeqAA.L;
        outScorePerCol = (diagLen > 0) ? ((float)idCnt / (float)diagLen) : 0.0f;
        if (par.addBacktrace) {
            backtrace.append(diagLen, 'M');
        }
        return Matcher::result_t(tSeqAA.getDbKey(), idCnt, queryCov, targetCov, seqId,
                                  std::numeric_limits<double>::max(), diagLen,
                                  qOff, qOff + diagLen - 1, qSeqAA.L,
                                  tOff, tOff + diagLen - 1, tSeqAA.L, backtrace);
    }

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
    outScorePerCol = (res.diagonalLen > 0) ? ((float)res.score / (float)res.diagonalLen) : 0.0f;

    double evalue = computeEvalue(res.score, hitsPerQuery, m_ln, c_ln, pfp);

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

    // MATCHA substitution matrix — resolve from bundled matrices if available
    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "MATCHA substitution matrix (--matcha) is required\n";
        EXIT(EXIT_FAILURE);
    }
    std::string teaMatData;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == par.teaMatrixFile) {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData,
                                   par.substitutionMatrices[i].subMatDataLen);
            char *serialized = BaseMatrix::serialize(par.substitutionMatrices[i].name, matrixData);
            teaMatData.assign(serialized);
            free(serialized);
            break;
        }
    }
    const char *teaMatSource = teaMatData.empty() ? par.teaMatrixFile.c_str() : teaMatData.c_str();
    SubstitutionMatrix subMatTea(teaMatSource, 1.0, par.scoreBias);

    // AA substitution matrix (from --sub-mat, weighted by --aa-weight)
    float aaFactor = par.teaWeight;
    SubstitutionMatrix subMatAA(par.scoringMatrixFile.values.aminoacid().c_str(), aaFactor, par.scoreBias);

    Debug::Progress progress(resultReader.getSize());

    double hitsPerQuery = static_cast<double>(par.maxResListLen);

    // E-value parameters (configurable via --loglinear-m, --loglinear-c, --p-fp)
    const double loglinearM_ln = par.loglinearM * 2.302585093;
    const double loglinearC_ln = par.loglinearC * 2.302585093;
    const double pfp = par.pFP;

    // --filter-hits + --rescore-mode: mirror mmseqs's rescorediagonal so the
    // linclust pipeline's per-step gating works as documented.
    float scorePerColThr = 0.0f;
    if (par.filterHits) {
        if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
            Debug(Debug::WARNING) << "HAMMING distance cannot be used to filter hits; switching to --rescore-mode 1\n";
            par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
        }
        std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
            ? std::string((const char*)CovSeqidQscPercMinDiag_lib, CovSeqidQscPercMinDiag_lib_len)
            : std::string((const char*)CovSeqidQscPercMinDiagTargetCov_lib, CovSeqidQscPercMinDiagTargetCov_lib_len);
        scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    }

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

                    float currScorePerCol = 0.0f;
                    Matcher::result_t res = ungappedAlignTea(qSeqAA, qSeqTea,
                                                              tSeqAA, tSeqTea,
                                                              static_cast<short>(prefHit.diagonal),
                                                              subMatAA, subMatTea, hitsPerQuery,
                                                              loglinearM_ln, loglinearC_ln, pfp,
                                                              backtrace, par, currScorePerCol);

                    if (res.dbKey == UINT_MAX) {
                        rejected++;
                        continue;
                    }

                    // mmseqs's --filter-hits gate: accept iff score-per-col passes
                    // the precision-library threshold, OR all of the existing
                    // (covThr / seqIdThr / alnLenThr / evalThr) criteria pass.
                    bool passedFilterHits = (par.filterHits && currScorePerCol >= scorePerColThr);
                    if (passedFilterHits ||
                        Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
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
