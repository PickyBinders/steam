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
#include "QueryMatcher.h"
#include "DatabaseDiversity.h"
#include "AlignmentSeedEvidence.h"
#include "seed_markov.txt.h"

#include <cmath>
#include <cstdlib>

#ifdef OPENMP
#include <omp.h>
#endif

// Continuous piecewise-loglinear E-value model.
// E(s) = P(FP) * searchSpaceSize * 10^(c + m_low*s
//                         + (m_high-m_low)*max(0,s-breakpoint)).
// The denominator is fixed MinHash diversity of the complete target database.
static double computeEvalue(double rawScore, double searchSpaceSize,
                            double mLow_ln, double mHigh_ln,
                            double breakpoint, double c_ln, double pfp) {
    double exponent = mLow_ln * rawScore + c_ln;
    if (breakpoint > 0.0 && rawScore > breakpoint) {
        exponent += (mHigh_ln - mLow_ln) * (rawScore - breakpoint);
    }
    return pfp * searchSpaceSize * exp(exponent);
}

static int doTeaAlign(TeaSmithWaterman &teaSW,
                      Sequence &tSeqAA, Sequence &tSeqTea,
                      unsigned int querySeqLen, unsigned int targetSeqLen,
                      double searchSpaceSize, double mLow_ln, double mHigh_ln,
                      double breakpoint, double c_ln, double pfp,
                      bool calibrated,
                      const AlignmentSeedEvidence *seedModel,
                      const AlignmentSeedEvidence::Query &seedQuery, const char *targetTea,
                      Matcher::result_t &res, std::string &backtrace,
                      const LocalParameters &par) {
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

    const int rankingScore = static_cast<int>(align.score1);
    align.evalue = (seedModel || !calibrated) ? std::numeric_limits<double>::infinity()
        : computeEvalue(rankingScore, searchSpaceSize, mLow_ln, mHigh_ln, breakpoint, c_ln, pfp);
    // Seed evidence is available only after traceback. Never reject using
    // the uncorrected score or the ranking-only infinity placeholder.
    if (calibrated && !seedModel && align.evalue > par.evalThr) {
        return -1;
    }

    // Full alignment with backtrace
    bool blockAlignFailed = false;
    if (teaSW.isProfileSearch() == false) {
        TeaSmithWaterman::s_align alignTmp = teaSW.alignStartPosBacktraceBlock(
            tSeqAA.numSequence, tSeqTea.numSequence, targetSeqLen,
            par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
            backtrace, align);
        // The block aligner returns its failure sentinel in a copy; preserve align for the fallback.
        if (alignTmp.score1 == UINT32_MAX) {
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

    int bitScore = rankingScore;
    if (seedModel) {
        bitScore = AlignmentSeedEvidence::score(rankingScore,
            seedModel->evidence(seedQuery, targetTea, targetSeqLen,
                                align.qStartPos1, align.dbStartPos1, backtrace));
    }
    align.evalue = !calibrated ? std::numeric_limits<double>::infinity()
        : computeEvalue(bitScore, searchSpaceSize, mLow_ln, mHigh_ln, breakpoint, c_ln, pfp);
    if (align.evalue > par.evalThr) {
        return -1;
    }

    res = Matcher::result_t(tSeqAA.getDbKey(), bitScore, align.qCov, align.tCov, seqId, align.evalue,
                            alnLength, align.qStartPos1, align.qEndPos1, querySeqLen,
                            align.dbStartPos1, align.dbEndPos1, targetSeqLen, backtrace);
    return 0;
}

int teaalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool calibrated = par.validateCalibration();
    const bool seedEnabled = par.seedCorrection;
    AlignmentSeedEvidence seedModel;
    if (seedEnabled) {
        if (par.compBiasCorrectionScale != 0.0f) {
            Debug(Debug::ERROR) << "Seed correction requires --comp-bias-corr-scale 0\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.seedModelFile == "seed_markov.txt") {
            seedModel.loadText(std::string((const char *)seed_markov_txt, seed_markov_txt_len), par.seedPattern);
        } else {
            seedModel.load(par.seedModelFile, par.seedPattern);
        }
        Debug(Debug::INFO) << "Alignment seed correction: floor(S*(1+e)+0.5)\n";
    }

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    bool sameDB = (par.db1.compare(par.db2) == 0);

    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;

    // Load TEA target database (main DB = TEA sequences)
    IndexReader tTeaDbr(par.db2, par.threads,
                        alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                        (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    // Load AA target companion database (_aa suffix)
    const std::string targetBase = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    std::string tAaDbName = targetBase + "_aa";
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
        std::string qAaDbName = PrefilteringIndexReader::dbPathWithoutIndex(par.db1) + "_aa";
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
    SubstitutionMatrix subMatTea(teaMatSource, par.teaScale, par.scoreBias);

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

    double searchSpaceSize = 0.0;
    if (!calibrated) {
        // Ranking-only custom scoring has no significance estimate.
        searchSpaceSize = 0.;
    } else {
        const size_t targetCount = tTeaDbr.sequenceReader->getSize();
        DatabaseDiversityMetadata diversity;
        std::string diversityError;
        if (!DatabaseDiversity::read(targetBase, diversity, diversityError)
                || diversity.sequences != targetCount) {
            if (diversityError.empty()) {
                diversityError = "metadata/target sequence-count mismatch";
            }
            Debug(Debug::ERROR)
                << "Cannot compute diversity-adjusted E-values: "
                << diversityError << "\n";
            EXIT(EXIT_FAILURE);
        }
        searchSpaceSize = diversity.effectiveTargets;
        Debug(Debug::INFO)
            << "E-value search space: diversity-adjusted MinHash = "
            << searchSpaceSize << " (" << diversity.method << ")\n";
    }

    // E-value parameters are supplied in log10 units and evaluated in ln.
    const double loglinearM_ln = par.loglinearM * 2.302585093;
    const double loglinearMHigh_ln = par.loglinearMHigh * 2.302585093;
    const double loglinearC_ln = par.loglinearC * 2.302585093;
    const double pfp = par.pFP;
    Debug(Debug::INFO) << "E-value params: m_low=" << par.loglinearM
                       << " m_high=" << par.loglinearMHigh
                       << " breakpoint=" << par.loglinearBreakpoint
                       << " c=" << par.loglinearC << " P(FP)=" << pfp << "\n";

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
                AlignmentSeedEvidence::Query seedQuery;
                if (seedEnabled) seedQuery = seedModel.prepare(querySeqTea, querySeqLen);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    const hit_t prefilterHit = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = prefilterHit.seqId;
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
                                   searchSpaceSize, loglinearM_ln,
                                   loglinearMHigh_ln, par.loglinearBreakpoint,
                                   loglinearC_ln, pfp,
                                   calibrated,
                                   seedEnabled ? &seedModel : NULL, seedQuery, targetSeqTea,
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

            // Sort by E-value ascending (= raw score descending).
            std::vector<size_t> order(alignmentResult.size());
            for (size_t i = 0; i < order.size(); i++) order[i] = i;
            if (order.size() > 1) {
                const auto &ar = alignmentResult;
                SORT_SERIAL(order.begin(), order.end(),
                    [&ar](size_t ia, size_t ib) {
                        const Matcher::result_t &a = ar[ia];
                        const Matcher::result_t &b = ar[ib];
                        if (a.eval != b.eval) return a.eval < b.eval;
                        if (a.score != b.score) return a.score > b.score;
                        return a.dbKey < b.dbKey;
                    });
            }
            for (size_t i = 0; i < order.size(); i++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[order[i]], par.addBacktrace);
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
