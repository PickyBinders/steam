#include "LocalParameters.h"
#include "Debug.h"
#include "matcha.out.h"
#include "FileUtil.h"
#include "IndexReader.h"
#include "Util.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TEA_WEIGHT(PARAM_TEA_WEIGHT_ID, "--aa-weight", "AA weight",
                         "Weight for amino acid substitution score in combined scoring (0.0 = structural alphabet only)",
                         typeid(float), (void *) &teaWeight,
                         "^[0-9]*(\\.[0-9]+)?$",
                         MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_TEA_MAT(PARAM_TEA_MAT_ID, "--matcha", "MATCHA substitution matrix",
                      "Path to MATCHA substitution matrix file",
                      typeid(std::string), (void *) &teaMatrixFile,
                      "",
                      MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_TEA_SCALE(PARAM_TEA_SCALE_ID, "--tea-scale", "TEA score scale",
                        "Alphabet score units per bit (2 for half-bit TEA matrices)",
                        typeid(float), (void *) &teaScale,
                        "^[0-9]*(\\.[0-9]+)?$",
                        MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOGLINEAR_M(PARAM_LOGLINEAR_M_ID, "--loglinear-m", "Low-score log-linear slope",
                          "Slope below the optional breakpoint in the continuous log-linear E-value model",
                          typeid(float), (void *) &loglinearM,
                          "^-?[0-9]*(\\.[0-9]+)?$",
                          MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOGLINEAR_M_HIGH(PARAM_LOGLINEAR_M_HIGH_ID, "--loglinear-m-high", "High-score log-linear slope",
                               "Slope above --loglinear-breakpoint; ignored when the breakpoint is zero",
                               typeid(float), (void *) &loglinearMHigh,
                               "^-?[0-9]*(\\.[0-9]+)?$",
                               MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOGLINEAR_BREAKPOINT(PARAM_LOGLINEAR_BREAKPOINT_ID, "--loglinear-breakpoint", "Log-linear score breakpoint",
                                   "Raw-score hinge for the calibrated MinHash E-value model; 0 with -e inf selects ranking-only output",
                                   typeid(float), (void *) &loglinearBreakpoint,
                                   "^[0-9]*(\\.[0-9]+)?$",
                                   MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOGLINEAR_C(PARAM_LOGLINEAR_C_ID, "--loglinear-c", "Log-linear intercept",
                          "Intercept c for the continuous log-linear E-value model",
                          typeid(float), (void *) &loglinearC,
                          "^-?[0-9]*(\\.[0-9]+)?$",
                          MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_P_FP(PARAM_P_FP_ID, "--p-fp", "P(FP) prior",
                   "Prior probability that a reported hit is a false positive",
                   typeid(float), (void *) &pFP,
                   "^[0-9]*(\\.[0-9]+)?$",
                   MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_UNGAPPED_TEA_AA(PARAM_UNGAPPED_TEA_AA_ID, "--ungapped-tea-aa", "Combined ungapped scoring",
                              "Rank prefilter hits by combined TEA+AA ungapped score before applying --max-seqs; unused in exhaustive search",
                              typeid(bool), (void *) &ungappedTeaAa,
                              "^[0-1]{1}$",
                              MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SEED_CORRECTION(PARAM_SEED_CORRECTION_ID, "--seed-correction", "Alignment seed correction",
                              "Apply alignment-supported query seed evidence once after traceback (0 for legacy alphabet scoring)",
                              typeid(bool), (void *) &seedCorrection, "^[0-1]{1}$",
                              MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SEED_MODEL(PARAM_SEED_MODEL_ID, "--seed-model", "Seed Markov model",
                         "Seed probability model; seed_markov.txt uses the bundled selected model",
                         typeid(std::string), (void *) &seedModelFile, "",
                         MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SEED_PATTERN(PARAM_SEED_PATTERN_ID, "--seed-pattern", "Alignment seed pattern",
                           "Alignment evidence pattern; search selects the same k5 or k6 pattern as its prefilter",
                           typeid(std::string), (void *) &seedPattern, "^(1101101|1110111)$",
                           MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_TFIDENT(PARAM_MIN_TFIDENT_ID, "--min-tfident", "Minimum TEA identity",
                          "Minimum fraction of identical TEA letters across alignment columns, including gaps",
                          typeid(float), (void *) &minTfident, "^(0(\\.[0-9]+)?|1(\\.0+)?)$",
                          MMseqsParameter::COMMAND_ALIGN)
{
    // Defaults
    teaWeight = 3.0;
    teaScale = 2.0;
    teaMatrixFile = "matcha.out";
    loglinearM = -0.005633978159801299;
    loglinearMHigh = -0.0014312498238873948;
    loglinearBreakpoint = 586.0;
    loglinearC = 0.2839183138420103;
    pFP = 1.0;
    ungappedTeaAa = true;
    seedCorrection = true;
    seedModelFile = "seed_markov.txt";
    seedPattern = "1101101";
    minTfident = 0.0f;
    // Register matcha.out as a bundled substitution matrix
    substitutionMatrices.push_back({"matcha.out", matcha_out, matcha_out_len});
    compBiasCorrection = 1;
    compBiasCorrectionScale = 0.0;
    maskMode = 0;
    exactKmerMatching = 1;
    gapOpen = MultiParam<NuclAA<int>>(NuclAA<int>(24, 5));
    gapExtend = MultiParam<NuclAA<int>>(NuclAA<int>(2, 2));
    maxResListLen = 2000;
    evalThr = 0.01;
    kmerSize = 5;
    spacedKmer = 1;
    spacedKmerPattern = "1101101";

    // Cluster workflows create align/rescorediagonal command lines from these
    // base parameter lists, not from teaalign/tearescorediagonal below.
    align.push_back(&PARAM_MIN_TFIDENT);
    rescorediagonal.push_back(&PARAM_MIN_TFIDENT);
    clusterworkflow.push_back(&PARAM_MIN_TFIDENT);
    linclustworkflow.push_back(&PARAM_MIN_TFIDENT);
    easyclusterworkflow.push_back(&PARAM_MIN_TFIDENT);
    easylinclustworkflow.push_back(&PARAM_MIN_TFIDENT);

    // createteadb parameters
    createteadb.push_back(&PARAM_WRITE_LOOKUP);
    createteadb.push_back(&PARAM_ID_OFFSET);
    createteadb.push_back(&PARAM_THREADS);
    createteadb.push_back(&PARAM_V);

    computediversity.push_back(&PARAM_THREADS);
    computediversity.push_back(&PARAM_V);

    createteasubdb = combineList(createsubdb, {&PARAM_THREADS});

    // teaalign = align + TEA-specific params + E-value params
    teaalign = combineList(align, {&PARAM_TEA_WEIGHT, &PARAM_TEA_MAT, &PARAM_TEA_SCALE,
                                    &PARAM_LOGLINEAR_M, &PARAM_LOGLINEAR_M_HIGH,
                                    &PARAM_LOGLINEAR_BREAKPOINT, &PARAM_LOGLINEAR_C, &PARAM_P_FP,
                                    &PARAM_SEED_CORRECTION, &PARAM_SEED_MODEL, &PARAM_SEED_PATTERN});

    // tearescorediagonal = mmseqs rescorediagonal + align + TEA-specific params.
    tearescorediagonal = combineList(rescorediagonal, align);
    tearescorediagonal = combineList(tearescorediagonal, {
        &PARAM_TEA_WEIGHT, &PARAM_TEA_MAT, &PARAM_TEA_SCALE,
        &PARAM_LOGLINEAR_M, &PARAM_LOGLINEAR_M_HIGH,
        &PARAM_LOGLINEAR_BREAKPOINT, &PARAM_LOGLINEAR_C, &PARAM_P_FP,
    });

    teaprefilter = combineList(prefilter, {&PARAM_TEA_MAT, &PARAM_TEA_WEIGHT,
        &PARAM_TEA_SCALE, &PARAM_UNGAPPED_TEA_AA});

    // teasearch = prefilter + teaalign + tearescorediagonal + common
    teasearchworkflow = combineList(prefilter, teaalign);
    teasearchworkflow = combineList(teasearchworkflow, tearescorediagonal);
    teasearchworkflow = combineList(teasearchworkflow, {&PARAM_UNGAPPED_TEA_AA,
                                                        &PARAM_RUNNER, &PARAM_REUSELATEST,
                                                        &PARAM_EXHAUSTIVE_SEARCH});

    // easyteasearch = teasearch + createteadb + convertalis
    easyteasearchworkflow = combineList(teasearchworkflow, createteadb);
    easyteasearchworkflow = combineList(easyteasearchworkflow, convertalignments);

    // Upstream clustering serializes these base lists for its child commands.
    align = teaalign;
    rescorediagonal = tearescorediagonal;
    clusterworkflow = combineList(clusterworkflow, teaalign);
    clusterworkflow = combineList(clusterworkflow, {&PARAM_UNGAPPED_TEA_AA});
    linclustworkflow = combineList(linclustworkflow, tearescorediagonal);
    linclustworkflow = combineList(linclustworkflow, {&PARAM_UNGAPPED_TEA_AA,
        &PARAM_SEED_CORRECTION, &PARAM_SEED_MODEL, &PARAM_SEED_PATTERN});
    easyclusterworkflow = combineList(easyclusterworkflow, clusterworkflow);
    easylinclustworkflow = combineList(easylinclustworkflow, linclustworkflow);
}

constexpr const char *STEAM_K5_PATTERN = "1101101";
constexpr const char *STEAM_K6_PATTERN = "1110111";

void LocalParameters::resolveKmerPattern(const std::string &targetDb) {
    LocalParameters &par = *this;
    if (par.PARAM_SPACED_KMER_PATTERN.wasSet && par.kmerSize == 0) {
        par.kmerSize = 0;
        for (char position : par.spacedKmerPattern) {
            par.kmerSize += position == '1';
        }
    }
    // Read only target index metadata, never target sequence data. Scale the
    // MMseqs k6 limit by one 20-letter seed position for STEAM's k5/k6 range.
    if (par.kmerSize == 0) {
        IndexReader target(targetDb, 1, IndexReader::SEQUENCES,
                           IndexReader::PRELOAD_NO, DBReader<unsigned int>::USE_INDEX);
        if (target.index != NULL) {
            par.kmerSize = PrefilteringIndexReader::getMetadata(target.index).kmerSize;
            if (!par.PARAM_SPACED_KMER_PATTERN.wasSet) {
                par.spacedKmerPattern = PrefilteringIndexReader::getSpacedPattern(target.index);
                par.spacedKmer = PrefilteringIndexReader::getMetadata(target.index).spacedKmer;
            }
        } else {
            const size_t k5Limit = IndexTable::getUpperBoundAACountForKmerSize(6) / 20;
            par.kmerSize = target.sequenceReader->getAminoAcidDBSize() < k5Limit ? 5 : 6;
        }
    }
    if (par.spacedKmer && !par.PARAM_SPACED_KMER_PATTERN.wasSet
            && par.spacedKmerPattern.empty()) {
        if (par.kmerSize == 5) {
            par.spacedKmerPattern = STEAM_K5_PATTERN;
        } else if (par.kmerSize == 6) {
            par.spacedKmerPattern = STEAM_K6_PATTERN;
        }
    }
    if (!par.spacedKmerPattern.empty()) {
        int patternWeight = 0;
        for (char position : par.spacedKmerPattern) {
            patternWeight += position == '1';
        }
        if (patternWeight != par.kmerSize) {
            Debug(Debug::ERROR)
                << "STEAM spaced k-mer pattern weight " << patternWeight
                << " does not match k=" << par.kmerSize << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO)
        << "STEAM k-mer selection: k=" << par.kmerSize
        << ", pattern=" << (par.spacedKmerPattern.empty()
                              ? "<contiguous>" : par.spacedKmerPattern)
        << "\n";
}


void LocalParameters::configureUngappedScoring(bool combined) {
    if (combined && teaScale != 2.0f) {
        Debug(Debug::ERROR) << "Combined ungapped scoring requires --tea-scale 2; use --ungapped-tea-aa 0 for other scales\n";
        EXIT(EXIT_FAILURE);
    }
    if (combined && !diagonalScoring) {
        Debug(Debug::ERROR) << "Combined ungapped scoring requires --diag-score 1\n";
        EXIT(EXIT_FAILURE);
    }
    setenv("STEAM_PREFILTER_CONFIGURED", "1", 1);
    setenv("STEAM_SPLIT_INVARIANT", "1", 1);
    setenv("STEAM_UNGAPPED_AA", combined ? "1" : "0", 1);
    setenv("STEAM_UNGAPPED_AA_SUBMAT", scoringMatrixFile.values.aminoacid().c_str(), 1);
    setenv("STEAM_UNGAPPED_AA_WEIGHT", SSTR(teaWeight).c_str(), 1);
}

void LocalParameters::configurePrefilter() {
    // Workflow children already receive a TEA --sub-mat and the separate AA
    // matrix. Do not reinterpret that TEA matrix as an AA scoring matrix.
    const char *configured = getenv("STEAM_PREFILTER_CONFIGURED");
    if (configured != NULL && std::string(configured) == "1") return;
    if (teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--matcha is required for STEAM prefilter\n";
        EXIT(EXIT_FAILURE);
    }
    for (const std::string &input : {db1, db2}) {
        const std::string base = PrefilteringIndexReader::dbPathWithoutIndex(input);
        if (!Parameters::isEqualDbtype(FileUtil::parseDbType(base.c_str()), DBTYPE_AMINO_ACIDS)
                || (ungappedTeaAa && !Parameters::isEqualDbtype(
                    FileUtil::parseDbType((base + "_aa").c_str()), DBTYPE_AMINO_ACIDS))) {
            Debug(Debug::ERROR) << "Combined STEAM prefilter requires paired TEA/AA sequence databases\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (!PARAM_K.wasSet) kmerSize = 0;
    if (!PARAM_SPACED_KMER_PATTERN.wasSet) spacedKmerPattern.clear();
    resolveKmerPattern(db2);
    configureUngappedScoring(ungappedTeaAa);
    scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(teaMatrixFile, teaMatrixFile));
    initMatrices();
    compBiasCorrectionScale = 0.0f;
}

void LocalParameters::configureClustering(bool linear) {
    if (!Parameters::isEqualDbtype(FileUtil::parseDbType(db1.c_str()), DBTYPE_AMINO_ACIDS)
            || !Parameters::isEqualDbtype(FileUtil::parseDbType((db1 + "_aa").c_str()), DBTYPE_AMINO_ACIDS)) {
        Debug(Debug::ERROR) << "STEAM clustering requires paired sequence databases\n";
        EXIT(EXIT_FAILURE);
    }
    if (teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--matcha is required for STEAM clustering\n";
        EXIT(EXIT_FAILURE);
    }
    // Linclust selects its own k for the reduced alphabet; it must not inherit
    // search's k5/k6 selection, even when the caller enables spaced k-mers.
    if (spacedKmer && !linear) {
        if (!PARAM_K.wasSet) kmerSize = 0;
        if (!PARAM_SPACED_KMER_PATTERN.wasSet) spacedKmerPattern.clear();
        resolveKmerPattern(db1);
        if (!PARAM_SEED_PATTERN.wasSet) seedPattern = spacedKmerPattern;
        if (seedCorrection && seedPattern != spacedKmerPattern) {
            Debug(Debug::ERROR) << "Alignment and prefilter seed patterns must agree\n";
            EXIT(EXIT_FAILURE);
        }
    }
    validateCalibration();
    setenv("MMSEQS_PREFILTER_SUBMAT", teaMatrixFile.c_str(), 1);
    configureUngappedScoring(ungappedTeaAa);
}

bool LocalParameters::validateCalibration() {
    const bool explicitFit = PARAM_LOGLINEAR_M.wasSet && PARAM_LOGLINEAR_M_HIGH.wasSet
        && PARAM_LOGLINEAR_C.wasSet && PARAM_LOGLINEAR_BREAKPOINT.wasSet;
    const std::string aaMatrix = scoringMatrixFile.values.aminoacid();
    const bool selected = seedCorrection && teaMatrixFile == "matcha.out"
        && seedModelFile == "seed_markov.txt" && teaWeight == 3.0f && teaScale == 2.0f
        && gapOpen.values.aminoacid() == 24 && gapExtend.values.aminoacid() == 2
        && compBiasCorrection == 1 && compBiasCorrectionScale == 0.0f && scoreBias == 0.0f
        && (aaMatrix == "blosum62.out" || aaMatrix.find("blosum62.out:") == 0);
    if (std::isnan(evalThr) || evalThr < 0) {
        Debug(Debug::ERROR) << "E-value threshold must be nonnegative\n";
        EXIT(EXIT_FAILURE);
    }
    if ((!selected && !explicitFit) || loglinearBreakpoint == 0.0f) {
        if (!std::isinf(evalThr)) {
            Debug(Debug::ERROR) << "Custom/legacy scoring requires -e inf for ranking-only search or all four --loglinear-* calibration parameters with a positive breakpoint\n";
            EXIT(EXIT_FAILURE);
        }
        loglinearBreakpoint = 0.0f;
        PARAM_LOGLINEAR_BREAKPOINT.wasSet = true;
        return false;
    }
    if (!std::isfinite(loglinearM) || !std::isfinite(loglinearMHigh)
            || !std::isfinite(loglinearC) || !std::isfinite(loglinearBreakpoint)
            || loglinearM >= 0 || loglinearMHigh >= 0 || loglinearBreakpoint <= 0) {
        Debug(Debug::ERROR) << "Calibration requires finite negative slopes and a positive breakpoint\n";
        EXIT(EXIT_FAILURE);
    }
    return true;
}
