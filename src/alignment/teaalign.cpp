#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "PairedClusterAlignment.h"
#include "BlockAligner.h"
#include <atomic>
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
#include <array>

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

using ResidueCounts = std::array<unsigned int, 256>;

static ResidueCounts countResidues(const unsigned char *sequence, unsigned int length) {
    ResidueCounts counts = {};
    for (unsigned int i = 0; i < length; ++i) ++counts[sequence[i]];
    return counts;
}

// Every identical aligned pair consumes one occurrence from each sequence.
static unsigned int sharedResidues(ResidueCounts remaining,
                                   const unsigned char *target, unsigned int length) {
    unsigned int shared = 0;
    for (unsigned int i = 0; i < length; ++i) {
        if (remaining[target[i]] > 0) {
            --remaining[target[i]];
            ++shared;
        }
    }
    return shared;
}

static unsigned int minimumCoveredLength(unsigned int length, float coverage) {
    unsigned int low = 0, high = length;
    // Match computeCov()'s float division, including threshold-rounding edges.
    while (low < high) {
        const unsigned int middle = low + (high - low) / 2;
        if (static_cast<float>(middle) / static_cast<float>(length) >= coverage)
            high = middle;
        else
            low = middle + 1;
    }
    return low;
}

static int doTeaAlign(TeaSmithWaterman &teaSW,
                      Sequence &tSeqAA, Sequence &tSeqTea,
                      unsigned int querySeqLen, unsigned int targetSeqLen,
                      double searchSpaceSize, double mLow_ln, double mHigh_ln,
                      double breakpoint, double c_ln, double pfp,
                      bool calibrated, bool isIdentity,
                      const AlignmentSeedEvidence *seedModel,
                      const AlignmentSeedEvidence::Query &seedQuery, double maxSeedEvidence,
                      const char *targetTea,
                      Matcher::result_t &res, std::string &backtrace,
                      const LocalParameters &par) {
    float seqId = 0.0;
    backtrace.clear();

    // Score + end position
    TeaSmithWaterman::s_align align = teaSW.alignScoreEndPos<TeaSmithWaterman::PROFILE>(
        tSeqAA.numSequence, tSeqTea.numSequence, targetSeqLen,
        par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
        querySeqLen / 2);

    if (teaSW.isProfileSearch() && (align.score1 == 0 || align.qEndPos1 < 0 || align.dbEndPos1 < 0)) {
        return -1;
    }

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
    if (calibrated && seedModel && std::isfinite(par.evalThr)) {
        // Skip traceback only if even supporting every query word cannot pass.
        // An unrepresentable hypothetical score must fall back to traceback,
        // not trigger score()'s guard when the actual score may still fit.
        const double upper = std::floor(static_cast<double>(rankingScore)
                                        * (1. + maxSeedEvidence) + .5);
        if (rankingScore >= 0 && std::isfinite(upper)
                && upper <= std::numeric_limits<int>::max()) {
            const int bestCorrected = AlignmentSeedEvidence::score(rankingScore, maxSeedEvidence);
            if (computeEvalue(bestCorrected, searchSpaceSize, mLow_ln, mHigh_ln,
                              breakpoint, c_ln, pfp) > par.evalThr) {
                return -1;
            }
        }
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
    if (par.minTfident > 0.0f &&
        (backtrace.empty() || static_cast<double>(align.identicalTeaCnt) <
            static_cast<double>(par.minTfident) * alnLength)) {
        return -1;
    }
    // These criteria were already authoritative in checkCriteria(). Apply
    // them before scanning/sorting alignment-supported seeds, preserving its
    // identity-hit exemption. The final check remains in place.
    if (!isIdentity && (!(seqId >= par.seqIdThr)
            || !Util::hasCoverage(par.covThr, par.covMode, align.qCov, align.tCov)
            || !Util::hasAlignmentLength(par.alnLenThr, alnLength))) {
        return -1;
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

static std::string teaMatrixSource(LocalParameters &par) {
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
    return teaMatData.empty() ? par.teaMatrixFile : teaMatData;
}

namespace {

class SteamClusterAlignment final : public PairedClusterAlignment {
public:
    explicit SteamClusterAlignment(LocalParameters &parameters)
        : par(parameters), calibrated(parameters.validateCalibration()),
          teaMatrix(teaMatrixSource(parameters).c_str(), par.teaScale, par.scoreBias),
          aaMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), par.teaWeight, par.scoreBias),
          tinyTea(flatten(teaMatrix)), tinyAA(flatten(aaMatrix)),
          queryAA(openAA(par.db1)),
          memberAA(par.filterSeqDBFile.empty() || par.filterSeqDBFile == par.db1
                   ? nullptr : openAA(par.filterSeqDBFile)),
          searchSpace(calibrated ? effectiveTargets(par.filterSeqDBFile.empty() ? par.db1 : par.filterSeqDBFile) : 0.0) {
        if (par.seedCorrection) {
            if (par.seedModelFile == "seed_markov.txt") {
                evidence.loadText(std::string(reinterpret_cast<const char *>(seed_markov_txt),
                                              seed_markov_txt_len), par.seedPattern);
            } else {
                evidence.load(par.seedModelFile, par.seedPattern);
            }
        }
        Debug(Debug::INFO) << "Linclust2 paired TEA+AA alignment, native banded extension, AA weight=" << par.teaWeight
                           << ", effective targets=" << searchSpace << "\n";
    }

    ~SteamClusterAlignment() override {
        Debug(Debug::INFO) << "Linclust2 paired alignment counts: ungapped accepted=" << ungappedAccepted.load()
                           << ", banded attempted=" << bandedTried.load() << ", banded accepted=" << bandedAccepted.load()
                           << ", exact fallback=" << exactTried.load() << "\n";
        if (memberAA) memberAA->close();
        queryAA->close();
    }

    std::unique_ptr<Worker> worker(unsigned int thread) override {
        return std::unique_ptr<Worker>(new Thread(*this, thread));
    }

private:
    static std::vector<int8_t> flatten(const SubstitutionMatrix &matrix) {
        std::vector<int8_t> values;
        for (int row = 0; row < matrix.alphabetSize; ++row) {
            values.insert(values.end(), matrix.subMatrix[row], matrix.subMatrix[row] + matrix.alphabetSize);
        }
        return values;
    }

    std::unique_ptr<DBReader<DBKeyType>> openAA(const std::string &base) const {
        const std::string name = base + "_aa";
        auto database = std::unique_ptr<DBReader<DBKeyType>>(new DBReader<DBKeyType>(
            name.c_str(), (name + ".index").c_str(), par.threads,
            DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX));
        database->open(DBReader<DBKeyType>::NOSORT);
        return database;
    }

    static double effectiveTargets(const std::string &base) {
        DatabaseDiversityMetadata diversity;
        std::string error;
        if (!DatabaseDiversity::read(base, diversity, error)) {
            Debug(Debug::ERROR) << "Cannot compute Linclust2 E-values: " << error << "\n";
            EXIT(EXIT_FAILURE);
        }
        return diversity.effectiveTargets;
    }

    class Thread final : public Worker {
    public:
        Thread(SteamClusterAlignment &context, unsigned int index)
            : ctx(context), thread(index),
              queryTea(ctx.par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &ctx.teaMatrix, 0, false, ctx.par.compBiasCorrection),
              queryAA(ctx.par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &ctx.aaMatrix, 0, false, ctx.par.compBiasCorrection),
              targetTea(ctx.par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &ctx.teaMatrix, 0, false, ctx.par.compBiasCorrection),
              targetAA(ctx.par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &ctx.aaMatrix, 0, false, ctx.par.compBiasCorrection) {
            const auto &par = ctx.par;
            const int open = par.gapOpen.values.aminoacid(), extend = par.gapExtend.values.aminoacid();
            if (open > extend && extend > 0 && open <= 127) {
                const size_t maxLength = ctx.memberAA ? std::max(ctx.queryAA->getMaxSeqLen(), ctx.memberAA->getMaxSeqLen())
                                                       : ctx.queryAA->getMaxSeqLen();
                banded.reset(new BlockAligner(Parameters::DBTYPE_AMINO_ACIDS, maxLength,
                    &ctx.aaMatrix, nullptr, nullptr, par.compBiasCorrection, par.compBiasCorrectionScale,
                    -open, -extend, &ctx.teaMatrix));
            }
        }

        ~Thread() override {
            ctx.ungappedAccepted.fetch_add(ungappedAccepted, std::memory_order_relaxed);
            ctx.bandedTried.fetch_add(bandedTried, std::memory_order_relaxed);
            ctx.bandedAccepted.fetch_add(bandedAccepted, std::memory_order_relaxed);
            ctx.exactTried.fetch_add(exactTried, std::memory_order_relaxed);
        }

        void prepare(DBKeyType key, const char *tea, unsigned int length) override {
            mapAA(*ctx.queryAA, queryAA, key, length);
            queryTea.mapSequence(0, key, tea, length);
            if (banded) banded->initQuery(&queryAA, &queryTea);
            exactReady = false;
            if (ctx.par.seedCorrection) seedQuery = ctx.evidence.prepare(tea, length);
            maxEvidence = ctx.par.seedCorrection && ctx.calibrated && std::isfinite(ctx.par.evalThr)
                          ? AlignmentSeedEvidence::maxEvidence(seedQuery) : 0.0;
        }

        bool align(DBKeyType key, const char *tea, unsigned int length,
                   unsigned short diagonal, bool member, Matcher::result_t &result) override {
            const auto &par = ctx.par;
            if (!Util::canBeCovered(par.covThr, par.covMode, queryTea.L, length)) return false;
            mapAA(ctx.memberAA ? *ctx.memberAA : *ctx.queryAA, targetAA, key, length);
            targetTea.mapSequence(0, key, tea, length);
            if (!banded) return exactAlign(tea, result);

            banded->setTargetCompanion(&targetTea);
            if (!member) memberDiagonal = diagonal;
            const auto seed = banded->pairedUngapped(&targetAA, memberDiagonal);
            const int qStart = seed.startPos + (seed.diagonal >= 0 ? seed.distToDiagonal : 0);
            const int tStart = seed.startPos + (seed.diagonal < 0 ? seed.distToDiagonal : 0);
            if (seed.score > 0) {
                const int span = seed.endPos - seed.startPos + 1;
                backtrace.assign(span, 'M');
                if (finish(seed.score, qStart, tStart, qStart + span - 1, tStart + span - 1, tea, result)) {
                    ++ungappedAccepted;
                    return true;
                }
            }
            // MMseqs2 gives original members an unrestricted alignment when the diagonal fails.
            if (member) return exactAlign(tea, result);
            if (seed.score == 0) return false;
            // Native Linclust2 anchors on three consecutive matches of the seed alphabet.
            for (int i = 1; i < seed.endPos - seed.startPos; ++i) {
                const auto *q = queryTea.numSequence + qStart + i;
                const auto *t = targetTea.numSequence + tStart + i;
                if (q[-1] != t[-1] || q[0] != t[0] || q[1] != t[1]) continue;
                backtrace.clear();
                int rawScore = 0;
                const int xDrop = 32 * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid();
                ++bandedTried;
                const auto alignment = banded->bandedalign(&targetAA, qStart + i, tStart + i,
                    backtrace, xDrop, par.covThr, par.covMode, &rawScore);
                if (alignment.evalue < 0 || backtrace.empty()) return false;
                if (finish(rawScore, alignment.qStartPos1, alignment.dbStartPos1,
                           alignment.qEndPos1, alignment.dbEndPos1, tea, result)) {
                    memberDiagonal = 0;
                    ++bandedAccepted;
                    return true;
                }
                return false;
            }
            return false;
        }

    private:
        bool finish(int rawScore, int qStart, int tStart, int qEnd, int tEnd,
                    const char *tea, Matcher::result_t &result) {
            const auto &par = ctx.par;
            const float qCov = SmithWaterman::computeCov(qStart, qEnd, queryAA.L);
            const float tCov = SmithWaterman::computeCov(tStart, tEnd, targetAA.L);
            const size_t span = backtrace.size();
            if (rawScore <= 0 || !Util::hasCoverage(par.covThr, par.covMode, qCov, tCov)
                    || !Util::hasAlignmentLength(par.alnLenThr, span)) return false;
            if (ctx.calibrated && std::isfinite(par.evalThr)
                    && evalue(std::floor(rawScore * (1.0 + maxEvidence) + .5)) > par.evalThr) return false;
            unsigned int aaIds = 0, teaIds = 0;
            int q = qStart, t = tStart;
            for (char op : backtrace) {
                if (op == 'M') {
                    aaIds += queryAA.numSequence[q] == targetAA.numSequence[t];
                    teaIds += queryTea.numSequence[q] == targetTea.numSequence[t];
                }
                q += op != 'D';
                t += op != 'I';
            }
            const float seqId = Util::computeSeqId(par.seqIdMode, aaIds, queryAA.L, targetAA.L, span);
            if (seqId < par.seqIdThr || teaIds < static_cast<double>(par.minTfident) * span) return false;
            const int score = par.seedCorrection ? AlignmentSeedEvidence::score(rawScore,
                ctx.evidence.evidence(seedQuery, tea, targetTea.L, qStart, tStart, backtrace)) : rawScore;
            const double eval = ctx.calibrated ? evalue(score) : std::numeric_limits<double>::infinity();
            if (eval > par.evalThr) return false;
            result = Matcher::result_t(targetAA.getDbKey(), score, qCov, tCov, seqId, eval, span,
                qStart, qEnd, queryAA.L, tStart, tEnd, targetAA.L, backtrace);
            return true;
        }

        double evalue(double score) const {
            const auto &par = ctx.par;
            return computeEvalue(score, ctx.searchSpace, par.loglinearM * lnTen, par.loglinearMHigh * lnTen,
                                 par.loglinearBreakpoint, par.loglinearC * lnTen, par.pFP);
        }

        bool exactAlign(const char *tea, Matcher::result_t &result) {
            ++exactTried;
            const auto &par = ctx.par;
            if (!exact) exact.reset(new TeaSmithWaterman(par.maxSeqLen, ctx.teaMatrix.alphabetSize,
                par.compBiasCorrection, par.compBiasCorrectionScale, &ctx.aaMatrix, &ctx.teaMatrix));
            if (!exactReady) {
                exact->ssw_init(&queryAA, &queryTea, ctx.tinyAA.data(), ctx.tinyTea.data(), &ctx.aaMatrix);
                exactReady = true;
            }
            return doTeaAlign(*exact, targetAA, targetTea, queryTea.L, targetTea.L, ctx.searchSpace,
                par.loglinearM * lnTen, par.loglinearMHigh * lnTen, par.loglinearBreakpoint,
                par.loglinearC * lnTen, par.pFP, ctx.calibrated, false,
                par.seedCorrection ? &ctx.evidence : nullptr, seedQuery, maxEvidence,
                tea, result, backtrace, par) == 0;
        }

        void mapAA(DBReader<DBKeyType> &database, Sequence &sequence, DBKeyType key,
                   unsigned int length) const {
            const auto id = database.getId(key);
            if (id == DB_ENTRY_NOT_FOUND || database.getSeqLen(id) != length) {
                Debug(Debug::ERROR) << "Linclust2 TEA/AA sequence key or length mismatch for " << key << "\n";
                EXIT(EXIT_FAILURE);
            }
            sequence.mapSequence(id, key, database.getData(id, thread), length);
        }

        SteamClusterAlignment &ctx;
        unsigned int thread;
        Sequence queryTea, queryAA, targetTea, targetAA;
        std::unique_ptr<BlockAligner> banded;
        std::unique_ptr<TeaSmithWaterman> exact;
        bool exactReady = false;
        unsigned short memberDiagonal = 0;
        double const lnTen = 2.302585092994046;
        size_t ungappedAccepted = 0, bandedTried = 0, bandedAccepted = 0, exactTried = 0;
        AlignmentSeedEvidence::Query seedQuery;
        double maxEvidence = 0.0;
        std::string backtrace;
    };

    const LocalParameters &par;
    bool calibrated;
    SubstitutionMatrix teaMatrix, aaMatrix;
    std::vector<int8_t> tinyTea, tinyAA;
    std::unique_ptr<DBReader<DBKeyType>> queryAA, memberAA;
    double searchSpace;
    AlignmentSeedEvidence evidence;
    std::atomic<size_t> ungappedAccepted{0}, bandedTried{0}, bandedAccepted{0}, exactTried{0};
};

} // namespace

PairedClusterAlignment *LocalParameters::createPairedClusterAlignment() {
    return new SteamClusterAlignment(*this);
}

int teaalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool profileQuery = Parameters::isEqualDbtype(FileUtil::parseDbType(
        PrefilteringIndexReader::dbPathWithoutIndex(par.db1).c_str()), Parameters::DBTYPE_HMM_PROFILE);
    if (profileQuery) {
        if (par.PARAM_TEA_WEIGHT.wasSet || par.PARAM_TEA_SCALE.wasSet ||
            (par.PARAM_SEED_CORRECTION.wasSet && par.seedCorrection) ||
            (par.PARAM_E.wasSet && par.evalThr != std::numeric_limits<double>::infinity()) ||
            par.seqIdThr != 0.0f || par.minTfident != 0.0f) {
            Debug(Debug::ERROR) << "Paired profile alignment uses preweighted rows and uncalibrated scores; "
                               << "sequence weighting, seed correction, E-value and identity filters do not apply\n";
            EXIT(EXIT_FAILURE);
        }
        par.evalThr = std::numeric_limits<double>::infinity();
        Debug(Debug::INFO) << "Paired profile alignment: raw ranking scores, no sequence calibration or seed correction\n";
    }
    const bool calibrated = !profileQuery && par.validateCalibration();
    const bool seedEnabled = !profileQuery && par.seedCorrection;
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

    const int queryTeaType = qTeaDbr->sequenceReader->getDbtype();
    const int queryAaType = qAaDbr->getDbtype();
    const int expectedQueryType = profileQuery ? Parameters::DBTYPE_HMM_PROFILE : Parameters::DBTYPE_AMINO_ACIDS;
    if (!Parameters::isEqualDbtype(queryTeaType, expectedQueryType) ||
        !Parameters::isEqualDbtype(queryAaType, expectedQueryType) ||
        !Parameters::isEqualDbtype(tTeaDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) ||
        !Parameters::isEqualDbtype(tAaDbr.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) ||
        qAaDbr->getSize() != qTeaDbr->sequenceReader->getSize() ||
        tAaDbr.getSize() != tTeaDbr.sequenceReader->getSize()) {
        Debug(Debug::ERROR) << "STEAM alignment requires matching paired queries and sequence targets\n";
        EXIT(EXIT_FAILURE);
    }
    if (profileQuery) {
        par.maxSeqLen = std::max<size_t>(par.maxSeqLen, qTeaDbr->sequenceReader->getMaxSeqLen());
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
    SubstitutionMatrix subMatTea(teaMatrixSource(par).c_str(), par.teaScale, par.scoreBias);

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

        Sequence qSeqAA(par.maxSeqLen, queryAaType, (const BaseMatrix *)&subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeqTea(par.maxSeqLen, queryTeaType, (const BaseMatrix *)&subMatTea, 0, false, par.compBiasCorrection);
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

                if (qAaDbr->getSeqLen(qAaDbr->getId(queryKey)) != querySeqLen) {
                    Debug(Debug::ERROR) << "Paired query sequence lengths differ\n";
                    EXIT(EXIT_FAILURE);
                }
                qSeqTea.mapSequence(id, queryKey, querySeqTea, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);

                // Unfiltered searches do no histogram work.
                const bool coverageColumns = par.covThr > 0. &&
                    (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL
                        || par.covMode == Parameters::COV_MODE_QUERY || par.covMode == Parameters::COV_MODE_TARGET);
                const bool columnBound = par.alnLenThr > 1 || coverageColumns;
                const bool aaIdentityBound = par.seqIdThr > 0.
                    && (columnBound || par.seqIdMode != Parameters::SEQ_ID_ALN_LEN);
                const bool teaIdentityBound = par.minTfident > 0. && columnBound;
                const bool identityBounds = aaIdentityBound || teaIdentityBound;
                ResidueCounts queryAACounts, queryTeaCounts;
                if (aaIdentityBound)
                    queryAACounts = countResidues(qSeqAA.numSequence, querySeqLen);
                if (teaIdentityBound)
                    queryTeaCounts = countResidues(qSeqTea.numSequence, querySeqLen);
                const unsigned int minimumQuerySpan = identityBounds &&
                    (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL || par.covMode == Parameters::COV_MODE_QUERY)
                    ? minimumCoveredLength(querySeqLen, par.covThr) : 0;

                teaSW.ssw_init(&qSeqAA, &qSeqTea, tinySubMatAA, tinySubMatTea, &subMatAA);
                AlignmentSeedEvidence::Query seedQuery;
                if (seedEnabled) seedQuery = seedModel.prepare(querySeqTea, querySeqLen);

                const double maxSeedEvidence = seedEnabled && calibrated && std::isfinite(par.evalThr)
                    ? AlignmentSeedEvidence::maxEvidence(seedQuery) : 0.;
                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    // Linclust's ungapped stage passes alignment rows here;
                    // ordinary search passes three-column prefilter rows.
                    // Both formats start with the target database key.
                    Util::parseKey(data, buffer);
                    const unsigned int dbKey = static_cast<unsigned int>(strtoul(buffer, NULL, 10));
                    data = Util::skipLine(data);
                    unsigned int targetId = tTeaDbr.sequenceReader->getId(dbKey);
                    const bool isIdentity = !profileQuery && (queryId == targetId && (par.includeIdentity || sameDB));

                    const int targetSeqLen = static_cast<int>(tTeaDbr.sequenceReader->getSeqLen(targetId));

                    if (tAaDbr.getSeqLen(tAaDbr.getId(dbKey)) != static_cast<size_t>(targetSeqLen)) {
                        Debug(Debug::ERROR) << "Paired target sequence lengths differ\n";
                        EXIT(EXIT_FAILURE);
                    }
                    if (Util::canBeCovered(par.covThr, par.covMode, static_cast<int>(querySeqLen), targetSeqLen) == false) {
                        rejected++;
                        continue;
                    }
                    char *targetSeqTea = tTeaDbr.sequenceReader->getData(targetId, thread_idx);
                    char *targetSeqAA = tAaDbr.getData(tAaDbr.getId(dbKey), thread_idx);
                    tSeqTea.mapSequence(targetId, dbKey, targetSeqTea, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);

                    if (identityBounds && !isIdentity) {
                        unsigned int minimumColumns = std::max(1, par.alnLenThr);
                        minimumColumns = std::max(minimumColumns, minimumQuerySpan);
                        if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL || par.covMode == Parameters::COV_MODE_TARGET)
                            minimumColumns = std::max(minimumColumns,
                                minimumCoveredLength(targetSeqLen, par.covThr));
                        if ((aaIdentityBound &&
                                !(Util::computeSeqId(par.seqIdMode,
                                    sharedResidues(queryAACounts, tSeqAA.numSequence, targetSeqLen),
                                    querySeqLen, targetSeqLen, minimumColumns) >= par.seqIdThr))
                                || (teaIdentityBound &&
                                    static_cast<double>(sharedResidues(queryTeaCounts, tSeqTea.numSequence, targetSeqLen))
                                        < static_cast<double>(par.minTfident) * minimumColumns)) {
                            rejected++;
                            continue;
                        }
                    }

                    Matcher::result_t res;
                    if (doTeaAlign(teaSW, tSeqAA, tSeqTea, querySeqLen, targetSeqLen,
                                   searchSpaceSize, loglinearM_ln,
                                   loglinearMHigh_ln, par.loglinearBreakpoint,
                                   loglinearC_ln, pfp,
                                   calibrated, isIdentity,
                                   seedEnabled ? &seedModel : NULL, seedQuery, maxSeedEvidence, targetSeqTea,
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
