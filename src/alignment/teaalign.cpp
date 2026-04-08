#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "TeaSmithWaterman.h"
#include "SubstitutionMatrix.h"
#include "FileUtil.h"
#include "FastSort.h"

#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

// Log-linear E-value model following Edgar & Sahakyan (2025).
// E(s) = P(FP) * (H/Q) * 10^(m*s + c)
// where s = raw alignment score, H/Q = total reported hits / number of queries.
// m, c, P(FP) are configurable via --loglinear-m, --loglinear-c, --p-fp.
static double computeEvalue(double rawScore, double hitsPerQuery, double m_ln, double c_ln, double pfp) {
    return pfp * hitsPerQuery * exp(m_ln * rawScore + c_ln);
}

static int doTeaAlign(TeaSmithWaterman &teaSW,
                      Sequence &tSeqAA, Sequence &tSeqTea,
                      unsigned int querySeqLen, unsigned int targetSeqLen,
                      double hitsPerQuery, double m_ln, double c_ln, double pfp,
                      Matcher::result_t &res, std::string &backtrace,
                      const Parameters &par) {
    float seqId = 0.0;
    backtrace.clear();

    // Score + end position
    TeaSmithWaterman::s_align align = teaSW.alignScoreEndPos<TeaSmithWaterman::PROFILE>(
        tSeqAA.numSequence, tSeqTea.numSequence, targetSeqLen,
        par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
        querySeqLen / 2);

    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, align.qCov, align.tCov));
    if (hasLowerCoverage) {
        return -1;
    }

    align.evalue = computeEvalue(align.score1, hitsPerQuery, m_ln, c_ln, pfp);
    if (align.evalue > par.evalThr) {
        return -1;
    }

    // Full alignment with backtrace
    bool blockAlignFailed = false;
    if (teaSW.isProfileSearch() == false) {
        TeaSmithWaterman::s_align alignTmp = teaSW.alignStartPosBacktraceBlock(
            tSeqAA.numSequence, tSeqTea.numSequence, targetSeqLen,
            par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
            backtrace, align);
        if (align.score1 == UINT32_MAX) {
            blockAlignFailed = true;
        } else {
            align = alignTmp;
        }
    }

    if (blockAlignFailed || teaSW.isProfileSearch()) {
        align = teaSW.alignStartPosBacktrace<TeaSmithWaterman::PROFILE>(
            tSeqAA.numSequence, tSeqTea.numSequence, targetSeqLen,
            par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
            par.alignmentMode, backtrace, align,
            par.covMode, par.covThr, querySeqLen / 2);
    }

    unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1,
                                                        align.dbStartPos1, align.dbEndPos1);
    if (backtrace.size() > 0) {
        alnLength = backtrace.size();
        seqId = Util::computeSeqId(par.seqIdMode, align.identicalAACnt, querySeqLen, targetSeqLen, alnLength);
    }

    int bitScore = static_cast<int>(align.score1);
    align.evalue = computeEvalue(align.score1, hitsPerQuery, m_ln, c_ln, pfp);

    res = Matcher::result_t(tSeqAA.getDbKey(), bitScore, align.qCov, align.tCov, seqId, align.evalue,
                            alnLength, align.qStartPos1, align.qEndPos1, querySeqLen,
                            align.dbStartPos1, align.dbEndPos1, targetSeqLen, backtrace);
    return 0;
}

int teaalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    bool sameDB = (par.db1.compare(par.db2) == 0);

    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;

    // Load TEA target database (main DB = TEA sequences)
    IndexReader tTeaDbr(par.db2, par.threads,
                        alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                        (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    // Load AA target companion database (_aa suffix)
    std::string tAaDbName = par.db2 + "_aa";
    DBReader<unsigned int> tAaDbr(tAaDbName.c_str(), (tAaDbName + ".index").c_str(), par.threads,
                                   DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tAaDbr.open(DBReader<unsigned int>::NOSORT);

    // Load query databases
    IndexReader *qTeaDbr = NULL;
    DBReader<unsigned int> *qAaDbr = NULL;
    if (sameDB) {
        qTeaDbr = &tTeaDbr;
        qAaDbr = &tAaDbr;
    } else {
        qTeaDbr = new IndexReader(par.db1, par.threads,
                                   alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                   (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        std::string qAaDbName = par.db1 + "_aa";
        qAaDbr = new DBReader<unsigned int>(qAaDbName.c_str(), (qAaDbName + ".index").c_str(), par.threads,
                                             DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        qAaDbr->open(DBReader<unsigned int>::NOSORT);
    }

    // Load prefilter results
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                         DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype = Parameters::DBTYPE_ALIGNMENT_RES;
    if (alignmentIsExtended) {
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
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

    Debug(Debug::INFO) << "TEA alignment with AA weight=" << aaFactor << "\n";
    Debug::Progress progress(resultReader.getSize());

    // Build tiny substitution matrices for profile creation
    int8_t *tinySubMatAA = (int8_t *)mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t *tinySubMatTea = (int8_t *)mem_align(ALIGN_INT, subMatTea.alphabetSize * 32);

    for (int i = 0; i < subMatTea.alphabetSize; i++) {
        for (int j = 0; j < subMatTea.alphabetSize; j++) {
            tinySubMatTea[i * subMatTea.alphabetSize + j] = subMatTea.subMatrix[i][j];
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }

    // Compute actual H/Q (hits per query) from prefilter results for E-value computation.
    size_t totalPrefilterHits = 0;
    size_t totalQueries = 0;
    for (size_t id = 0; id < resultReader.getSize(); id++) {
        char *data = resultReader.getData(id, 0);
        if (*data != '\0') {
            totalQueries++;
            while (*data != '\0') {
                data = Util::skipLine(data);
                totalPrefilterHits++;
            }
        }
    }
    double hitsPerQuery = (totalQueries > 0) ? static_cast<double>(totalPrefilterHits) / static_cast<double>(totalQueries) : 500.0;
    Debug(Debug::INFO) << "H/Q (hits per query) = " << hitsPerQuery << " (" << totalPrefilterHits << " / " << totalQueries << ")\n";

    // E-value parameters (configurable via --loglinear-m, --loglinear-c, --p-fp)
    const double loglinearM_ln = par.loglinearM * 2.302585093;  // convert log10 to ln
    const double loglinearC_ln = par.loglinearC * 2.302585093;
    const double pfp = par.pFP;
    Debug(Debug::INFO) << "E-value params: m=" << par.loglinearM << " c=" << par.loglinearC << " P(FP)=" << pfp << "\n";

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        TeaSmithWaterman teaSW(par.maxSeqLen, subMatTea.alphabetSize, par.compBiasCorrection,
                                par.compBiasCorrectionScale, &subMatAA, &subMatTea);

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
                unsigned int queryId = qTeaDbr->sequenceReader->getId(queryKey);

                char *querySeqAA = qAaDbr->getData(qAaDbr->getId(queryKey), thread_idx);
                char *querySeqTea = qTeaDbr->sequenceReader->getData(queryId, thread_idx);
                unsigned int querySeqLen = qTeaDbr->sequenceReader->getSeqLen(queryId);

                qSeqTea.mapSequence(id, queryKey, querySeqTea, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);

                teaSW.ssw_init(&qSeqAA, &qSeqTea, tinySubMatAA, tinySubMatTea, &subMatAA);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int)strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tTeaDbr.sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB));

                    char *targetSeqTea = tTeaDbr.sequenceReader->getData(targetId, thread_idx);
                    char *targetSeqAA = tAaDbr.getData(tAaDbr.getId(dbKey), thread_idx);
                    const int targetSeqLen = static_cast<int>(tTeaDbr.sequenceReader->getSeqLen(targetId));

                    tSeqTea.mapSequence(targetId, dbKey, targetSeqTea, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);

                    if (Util::canBeCovered(par.covThr, par.covMode, static_cast<int>(querySeqLen), targetSeqLen) == false) {
                        rejected++;
                        continue;
                    }

                    Matcher::result_t res;
                    if (doTeaAlign(teaSW, tSeqAA, tSeqTea, querySeqLen, targetSeqLen,
                                   hitsPerQuery, loglinearM_ln, loglinearC_ln, pfp,
                                   res, backtrace, par) == -1) {
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
                // Sort by E-value ascending (= raw score descending)
                SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(),
                    [](const Matcher::result_t &a, const Matcher::result_t &b) {
                        if (a.eval != b.eval) return a.eval < b.eval;
                        if (a.score != b.score) return a.score > b.score;
                        return a.dbKey < b.dbKey;
                    });
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

    free(tinySubMatAA);
    free(tinySubMatTea);

    dbw.close();
    resultReader.close();
    tAaDbr.close();

    if (sameDB == false) {
        delete qTeaDbr;
        qAaDbr->close();
        delete qAaDbr;
    }

    return EXIT_SUCCESS;
}
