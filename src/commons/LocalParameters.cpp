#include "LocalParameters.h"
#include "Debug.h"
#include "matcha.out.h"
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
                        MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
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
                           MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT)
{
    // Defaults
    teaWeight = 3.0;
    teaScale = 2.0;
    teaMatrixFile = "matcha.out";
    loglinearM = -0.005628286904365784;
    loglinearMHigh = -0.0013916072449292318;
    loglinearBreakpoint = 578.0;
    loglinearC = 0.24858538629087792;
    pFP = 1.0;
    ungappedTeaAa = true;
    seedCorrection = true;
    seedModelFile = "seed_markov.txt";
    seedPattern = "1101101";
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

    // teaprefilter = base prefilter list (for the steam-specific prefilter alias)
    teaprefilter = prefilter;

    // teasearch = prefilter + teaalign + tearescorediagonal + common
    teasearchworkflow = combineList(prefilter, teaalign);
    teasearchworkflow = combineList(teasearchworkflow, tearescorediagonal);
    teasearchworkflow = combineList(teasearchworkflow, {&PARAM_UNGAPPED_TEA_AA,
                                                        &PARAM_RUNNER, &PARAM_REUSELATEST,
                                                        &PARAM_EXHAUSTIVE_SEARCH});

    // easyteasearch = teasearch + createteadb + convertalis
    easyteasearchworkflow = combineList(teasearchworkflow, createteadb);
    easyteasearchworkflow = combineList(easyteasearchworkflow, convertalignments);
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
