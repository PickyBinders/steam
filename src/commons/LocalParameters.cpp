#include "LocalParameters.h"

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TEA_WEIGHT(PARAM_TEA_WEIGHT_ID, "--tea-weight", "TEA AA weight",
                         "Weight for amino acid substitution score in TEA combined scoring (0.0 = TEA only)",
                         typeid(float), (void *) &teaWeight,
                         "^[0-9]*(\\.[0-9]+)?$",
                         MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_TEA_MAT(PARAM_TEA_MAT_ID, "--tea-mat", "TEA substitution matrix",
                      "Path to TEA substitution matrix file",
                      typeid(std::string), (void *) &teaMatrixFile,
                      "",
                      MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT)
{
    // Defaults
    teaWeight = 1.4;
    teaMatrixFile = "";

    // createteadb parameters
    createteadb.push_back(&PARAM_WRITE_LOOKUP);
    createteadb.push_back(&PARAM_ID_OFFSET);
    createteadb.push_back(&PARAM_THREADS);
    createteadb.push_back(&PARAM_V);

    // teaalign = align + TEA-specific params
    teaalign = combineList(align, {&PARAM_TEA_WEIGHT, &PARAM_TEA_MAT});

    // teasearch = prefilter + teaalign + common
    teasearchworkflow = combineList(prefilter, teaalign);
    teasearchworkflow = combineList(teasearchworkflow, {&PARAM_RUNNER, &PARAM_REUSELATEST});

    // easyteasearch = teasearch + createteadb + convertalis
    easyteasearchworkflow = combineList(teasearchworkflow, createteadb);
    easyteasearchworkflow = combineList(easyteasearchworkflow, convertalignments);

    // teacluster = prefilter + teaalign + rescorediagonal + clust + linclust + common
    teaclusterworkflow = combineList(prefilter, teaalign);
    teaclusterworkflow = combineList(teaclusterworkflow, rescorediagonal);
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
