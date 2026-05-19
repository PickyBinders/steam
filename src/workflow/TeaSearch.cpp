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
    cmd.addVariable("EXHAUSTIVE", par.exhaustiveSearch ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    // Override --sub-mat with --matcha for prefiltering (TEA k-mer matching)
    // Use lower comp bias scale for prefilter (matching foldseek)
    auto origScoringMatrixFile = par.scoringMatrixFile;
    par.scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(par.teaMatrixFile, par.teaMatrixFile));
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.teaprefilter).c_str());

    // Exhaustive search: mirror mmseqs2's Search.cpp logic.
    // teasearch.sh's fake_pref creates a synthetic prefilter that pairs every
    // query with every target, so the k-mer prefilter is bypassed entirely.
    // We still adjust covMode / maxResListLen / evalThr the same way mmseqs
    // does so alignment-stage filtering behaves correctly on the all-pairs
    // input.
    if (par.exhaustiveSearch) {
        const size_t queryDbSize  = FileUtil::countLines(par.db1Index.c_str());
        const size_t targetDbSize = FileUtil::countLines(par.db2Index.c_str());
        par.covMode = Util::swapCoverageMode(par.covMode);
        par.maxResListLen = std::max((size_t)300, queryDbSize);
        par.evalThr *= ((float) queryDbSize) / targetDbSize;
    }

    // Restore original --sub-mat so rescorediagonal and alignment use the AA matrix (e.g. BLOSUM62),
    // not the TEA matrix, for the amino acid scoring component
    par.scoringMatrixFile = origScoringMatrixFile;

    // Gapped alignment uses the user's E-value threshold
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.teaalign).c_str());

    std::string program(tmpDir + "/teasearch.sh");
    FileUtil::writeFile(program, teasearch_sh, teasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
