#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easyteacluster.sh.h"

int easyteacluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.teaclusterworkflow));
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

    // Input: tea_fasta aa_fasta
    cmd.addVariable("INPUT_TEA", par.filenames[0].c_str());
    cmd.addVariable("INPUT_AA", par.filenames[1].c_str());

    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.teaclusterworkflow, true).c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createteadb).c_str());
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    std::string program(tmpDir + "/easyteacluster.sh");
    FileUtil::writeFile(program, easyteacluster_sh, easyteacluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
