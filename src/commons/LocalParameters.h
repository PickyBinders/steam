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
    std::string teaMatrixFile;
    int minUngappedScore;

    PARAMETER(PARAM_TEA_WEIGHT)
    PARAMETER(PARAM_TEA_MAT)
    PARAMETER(PARAM_MIN_UNGAPPED_SCORE)

    // Parameter vectors for TEA commands
    std::vector<MMseqsParameter*> createteadb;
    std::vector<MMseqsParameter*> teaalign;
    std::vector<MMseqsParameter*> tearescorediagonal;

    // Workflow parameter vectors
    std::vector<MMseqsParameter*> teasearchworkflow;
    std::vector<MMseqsParameter*> easyteasearchworkflow;
    std::vector<MMseqsParameter*> teaclusterworkflow;
    std::vector<MMseqsParameter*> easyteaclusterworkflow;

    LocalParameters();

private:
    LocalParameters(const LocalParameters&);
    void operator=(const LocalParameters&);
};

#endif // STEAM_LOCAL_PARAMETERS_H
