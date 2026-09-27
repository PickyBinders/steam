#ifndef STEAM_LOCAL_PARAMETERS_H
#define STEAM_LOCAL_PARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initParameterSingleton();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    // TEA-specific output format codes
    static const int OUTFMT_TFIDENT = 100;   // TEA fractional identity
    static const int OUTFMT_TPIDENT = 101;   // TEA percent identity
    static const int OUTFMT_QTEASEQ = 102;   // query TEA full sequence
    static const int OUTFMT_TTEASEQ = 103;   // target TEA full sequence
    static const int OUTFMT_QTEAALN = 104;   // query TEA aligned sequence
    static const int OUTFMT_TTEAALN = 105;   // target TEA aligned sequence

    // TEA-specific parameters
    float teaWeight;
    float teaScale;
    std::string teaMatrixFile;
    float loglinearM;
    float loglinearMHigh;
    float loglinearBreakpoint;
    float loglinearC;
    float pFP;
    bool ungappedTeaAa;
    bool seedCorrection;
    std::string seedModelFile;
    std::string seedPattern;

    PARAMETER(PARAM_TEA_WEIGHT)
    PARAMETER(PARAM_TEA_MAT)
    PARAMETER(PARAM_TEA_SCALE)
    PARAMETER(PARAM_LOGLINEAR_M)
    PARAMETER(PARAM_LOGLINEAR_M_HIGH)
    PARAMETER(PARAM_LOGLINEAR_BREAKPOINT)
    PARAMETER(PARAM_LOGLINEAR_C)
    PARAMETER(PARAM_P_FP)
    PARAMETER(PARAM_UNGAPPED_TEA_AA)
    PARAMETER(PARAM_SEED_CORRECTION)
    PARAMETER(PARAM_SEED_MODEL)
    PARAMETER(PARAM_SEED_PATTERN)

    // Parameter vectors for TEA commands
    std::vector<MMseqsParameter*> createteadb;
    std::vector<MMseqsParameter*> computediversity;
    std::vector<MMseqsParameter*> createteasubdb;
    std::vector<MMseqsParameter*> teaalign;
    std::vector<MMseqsParameter*> tearescorediagonal;
    std::vector<MMseqsParameter*> teaprefilter;

    // Workflow parameter vectors
    std::vector<MMseqsParameter*> teasearchworkflow;
    std::vector<MMseqsParameter*> easyteasearchworkflow;

    LocalParameters();
    bool validateCalibration();

private:
    LocalParameters(const LocalParameters&);
    void operator=(const LocalParameters&);
};

#endif // STEAM_LOCAL_PARAMETERS_H
