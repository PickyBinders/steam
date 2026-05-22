// Steam wrappers for cluster / linclust / easy-cluster / easy-linclust.
//
// cluster + linclust: thin wrappers that setenv MMSEQS_PREFILTER_SUBMAT so
// the patched mmseqs Cluster.cpp/Linclust.cpp swap par.scoringMatrixFile to
// matcha.out around prefilter param string generation. The actual workflow
// runs upstream's clusteringworkflow() / linclust() entry points.
//
// easy-cluster + easy-linclust: drive our own easyteacluster.sh because
// steam's `createdb` takes <tea_fa> <aa_fa> <out>, not a variadic FASTA list,
// so mmseqs's easycluster.sh signature doesn't fit. The shell still delegates
// the heavy lifting to mmseqs's `cluster` / `linclust` (via our shadows).

#include <cassert>
#include <cstdlib>

#include "Command.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "easyteacluster.sh.h"

static void setupTeaEnv() {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setenv("MMSEQS_PREFILTER_SUBMAT", par.teaMatrixFile.empty() ? "matcha.out" : par.teaMatrixFile.c_str(), 1);
}

int teacluster(int argc, const char **argv, const Command &command) {
    setupTeaEnv();
    return clusteringworkflow(argc, argv, command);
}

int teaLinclust(int argc, const char **argv, const Command &command) {
    setupTeaEnv();
    return linclust(argc, argv, command);
}

static int runEasyTeaCluster(int argc, const char **argv, const Command &command,
                             const char *clusterModule,
                             const std::vector<MMseqsParameter*> &clusterWorkflowParams) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, clusterWorkflowParams));
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

    // Remaining filenames: <teaFasta> <aaFasta>
    if (par.filenames.size() < 2) {
        Debug(Debug::ERROR) << "easy-cluster requires two FASTA inputs: <teaFasta> <aaFasta>\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("INPUT_TEA", par.filenames[0].c_str());
    cmd.addVariable("INPUT_AA", par.filenames[1].c_str());

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createteadb).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(clusterWorkflowParams, true).c_str());
    cmd.addVariable("CLUSTER_MODULE", clusterModule);
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyteacluster.sh";
    FileUtil::writeFile(program, easyteacluster_sh, easyteacluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    assert(false);
    return EXIT_FAILURE;
}

int easyteacluster(int argc, const char **argv, const Command &command) {
    setupTeaEnv();
    LocalParameters &par = LocalParameters::getLocalInstance();
    return runEasyTeaCluster(argc, argv, command, "cluster", par.clusterworkflow);
}

int easytealinclust(int argc, const char **argv, const Command &command) {
    setupTeaEnv();
    LocalParameters &par = LocalParameters::getLocalInstance();
    return runEasyTeaCluster(argc, argv, command, "linclust", par.linclustworkflow);
}
