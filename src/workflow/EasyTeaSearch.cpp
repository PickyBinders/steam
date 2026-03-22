#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easyteasearch.sh.h"

static void teaSearchDefault(LocalParameters &par) {
    par.compBiasCorrectionScale = 0.5;
}

int easyteasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    teaSearchDefault(par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.teasearchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    // Query input: tea_fasta aa_fasta
    cmd.addVariable("QUERY_TEA", par.filenames[0].c_str());
    cmd.addVariable("QUERY_AA", par.filenames[1].c_str());

    // Target: can be tea_fasta+aa_fasta pair or pre-built DB
    // After popping tmpDir and RESULTS, filenames has: [query_tea, query_aa, target_tea, target_aa] or [query_tea, query_aa, targetDB]
    if (par.filenames.size() >= 4) {
        // Target is fasta pair: tea_fasta aa_fasta
        cmd.addVariable("TARGET_TEA", par.filenames[2].c_str());
        cmd.addVariable("TARGET_AA", par.filenames[3].c_str());
        cmd.addVariable("TARGET", par.filenames[2].c_str()); // will be overwritten by script
    } else {
        // Target is pre-built DB
        cmd.addVariable("TARGET", par.filenames[2].c_str());
    }

    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.teasearchworkflow, true).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createteadb).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    std::string program(tmpDir + "/easyteasearch.sh");
    FileUtil::writeFile(program, easyteasearch_sh, easyteasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
