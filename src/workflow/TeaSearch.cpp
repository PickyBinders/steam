#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "teasearch.sh.h"

static void teaSearchDefault(LocalParameters &par) {
    par.compBiasCorrectionScale = 0.5;
    par.maskMode = 0;
}

int teasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    teaSearchDefault(par);

    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--matcha is required for steam search\n";
        EXIT(EXIT_FAILURE);
    }

    std::string tmpDir = par.filenames.back();
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("QUERY", par.filenames[0].c_str());
    cmd.addVariable("TARGET", par.filenames[1].c_str());
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    // Override --sub-mat with --matcha for prefiltering (TEA k-mer matching)
    // Use lower comp bias scale for prefilter (matching foldseek)
    auto origScoringMatrixFile = par.scoringMatrixFile;
    par.scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(par.teaMatrixFile, par.teaMatrixFile));
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());

    // Restore original --sub-mat so rescorediagonal and alignment use the AA matrix (e.g. BLOSUM62),
    // not the TEA matrix, for the amino acid scoring component
    par.scoringMatrixFile = origScoringMatrixFile;

    // Rescorediagonal and alignment use --matcha (TEA) + --sub-mat (AA, with --aa-weight scaling)
    cmd.addVariable("RESCOREDIAGONAL_PAR", par.createParameterString(par.tearescorediagonal).c_str());
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.teaalign).c_str());

    std::string program(tmpDir + "/teasearch.sh");
    FileUtil::writeFile(program, teasearch_sh, teasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
