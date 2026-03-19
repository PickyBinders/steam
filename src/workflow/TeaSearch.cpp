#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "teasearch.sh.h"

int teasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--tea-mat is required for steam search\n";
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

    // Override --sub-mat with --tea-mat for prefiltering (TEA k-mer matching)
    par.scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(par.teaMatrixFile, par.teaMatrixFile));
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());

    // Alignment uses --tea-mat (TEA) + --sub-mat (AA, with --tea-weight scaling)
    // teaalign reads --tea-mat and --sub-mat directly from its own parameters
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.teaalign).c_str());

    std::string program(tmpDir + "/teasearch.sh");
    FileUtil::writeFile(program, teasearch_sh, teasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
