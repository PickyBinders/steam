#include "LocalParameters.h"
#include "matcha.out.h"

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
        PARAM_LOGLINEAR_M(PARAM_LOGLINEAR_M_ID, "--loglinear-m", "Log-linear slope",
                          "Slope parameter m for log-linear E-value model: E = P(FP) * (H/Q) * 10^(m*s+c)",
                          typeid(float), (void *) &loglinearM,
                          "^-?[0-9]*(\\.[0-9]+)?$",
                          MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOGLINEAR_C(PARAM_LOGLINEAR_C_ID, "--loglinear-c", "Log-linear intercept",
                          "Intercept parameter c for log-linear E-value model: E = P(FP) * (H/Q) * 10^(m*s+c)",
                          typeid(float), (void *) &loglinearC,
                          "^-?[0-9]*(\\.[0-9]+)?$",
                          MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_P_FP(PARAM_P_FP_ID, "--p-fp", "P(FP) prior",
                   "Prior probability that a reported hit is a false positive (0.5 = filtered, 1.0 = unfiltered)",
                   typeid(float), (void *) &pFP,
                   "^[0-9]*(\\.[0-9]+)?$",
                   MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT)
{
    // Defaults
    teaWeight = 1.4;
    teaMatrixFile = "matcha.out";
    loglinearM = -0.0182549591;
    loglinearC = 0.03214628;
    pFP = 1.0;
    // Register matcha.out as a bundled substitution matrix
    substitutionMatrices.push_back({"matcha.out", matcha_out, matcha_out_len});
    compBiasCorrection = 0;
    compBiasCorrectionScale = 0.5;
    maskMode = 0;
    exactKmerMatching = 1;
    gapOpen = MultiParam<NuclAA<int>>(NuclAA<int>(14, 5));
    gapExtend = MultiParam<NuclAA<int>>(NuclAA<int>(2, 2));
    maxResListLen = 2000;
    evalThr = 100.0;
    kmerSize = 6;
    spacedKmer = 1;
    spacedKmerPattern = "110101101";

    // createteadb parameters
    createteadb.push_back(&PARAM_WRITE_LOOKUP);
    createteadb.push_back(&PARAM_ID_OFFSET);
    createteadb.push_back(&PARAM_THREADS);
    createteadb.push_back(&PARAM_V);

    // teaalign = align + TEA-specific params + E-value params
    teaalign = combineList(align, {&PARAM_TEA_WEIGHT, &PARAM_TEA_MAT,
                                    &PARAM_LOGLINEAR_M, &PARAM_LOGLINEAR_C, &PARAM_P_FP});

    // tearescorediagonal = align + TEA-specific params (for ungapped diagonal rescoring)
    tearescorediagonal = combineList(align, {&PARAM_TEA_WEIGHT, &PARAM_TEA_MAT});

    // teasearch = prefilter + teaalign + tearescorediagonal + common
    teasearchworkflow = combineList(prefilter, teaalign);
    teasearchworkflow = combineList(teasearchworkflow, tearescorediagonal);
    teasearchworkflow = combineList(teasearchworkflow, {&PARAM_RUNNER, &PARAM_REUSELATEST, &PARAM_EXHAUSTIVE_SEARCH});

    // easyteasearch = teasearch + createteadb + convertalis
    easyteasearchworkflow = combineList(teasearchworkflow, createteadb);
    easyteasearchworkflow = combineList(easyteasearchworkflow, convertalignments);

    // teacluster = prefilter + teaalign + tearescorediagonal + clust + linclust + common
    teaclusterworkflow = combineList(prefilter, teaalign);
    teaclusterworkflow = combineList(teaclusterworkflow, tearescorediagonal);
    teaclusterworkflow = combineList(teaclusterworkflow, clust);
    teaclusterworkflow.push_back(&PARAM_CASCADED);
    teaclusterworkflow.push_back(&PARAM_CLUSTER_STEPS);
    teaclusterworkflow.push_back(&PARAM_CLUSTER_REASSIGN);
    teaclusterworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    teaclusterworkflow.push_back(&PARAM_REUSELATEST);
    teaclusterworkflow.push_back(&PARAM_RUNNER);
    teaclusterworkflow = combineList(teaclusterworkflow, linclustworkflow);

    // easyteacluster = teacluster + createteadb + result2repseq
    easyteaclusterworkflow = combineList(teaclusterworkflow, createteadb);
    easyteaclusterworkflow = combineList(easyteaclusterworkflow, result2repseq);
}
