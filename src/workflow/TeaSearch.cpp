#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "IndexReader.h"

#include "teasearch.sh.h"

static void teaSearchDefault(LocalParameters &par) {
        par.compBiasCorrectionScale = 0.0;
        par.maskMode = 0;
        par.kmerSize = 0;
        par.spacedKmerPattern = "";
    }


int teasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    teaSearchDefault(par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--matcha is required for steam search\n";
        EXIT(EXIT_FAILURE);
    }
    for (int db = 0; db < 2; ++db) {
        const std::string sequenceDb = PrefilteringIndexReader::dbPathWithoutIndex(par.filenames[db]);
        if (!Parameters::isEqualDbtype(FileUtil::parseDbType(sequenceDb.c_str()),
                                      Parameters::DBTYPE_AMINO_ACIDS)) {
            Debug(Debug::ERROR) << "STEAM search requires paired sequence databases\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.validateCalibration();
    par.resolveKmerPattern(par.filenames[1]);
    if (!par.PARAM_SEED_PATTERN.wasSet) par.seedPattern = par.spacedKmerPattern;
    if (par.seedCorrection && par.seedPattern != par.spacedKmerPattern) {
        Debug(Debug::ERROR) << "Alignment and prefilter seed patterns must agree\n";
        EXIT(EXIT_FAILURE);
    }
    // Exhaustive alignment bypasses candidate selection, not final score terms.
    const bool combinedUngapped = par.ungappedTeaAa && !par.exhaustiveSearch;
    par.configureUngappedScoring(combinedUngapped);

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
    // Exact seeds and combined AA+TEA diagonal scoring need no profile bias.
    auto origScoringMatrixFile = par.scoringMatrixFile;
    const float alignmentCompBiasCorrectionScale = par.compBiasCorrectionScale;
    par.scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(par.teaMatrixFile, par.teaMatrixFile));
    par.compBiasCorrectionScale = 0.0;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.teaprefilter).c_str());

    // Unlike MMseqs' transposed exhaustive workflow, the exhaustive path retain query
    // and target orientation. Keep coverage modes and the user's E cutoff.

    // Restore original --sub-mat so rescorediagonal and alignment use the AA matrix (e.g. BLOSUM62),
    // not the TEA matrix, for the amino acid scoring component
    par.scoringMatrixFile = origScoringMatrixFile;
    par.compBiasCorrectionScale = alignmentCompBiasCorrectionScale;

    // Gapped alignment uses the user's E-value threshold
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.teaalign).c_str());

    std::string program(tmpDir + "/teasearch.sh");
    FileUtil::writeFile(program, teasearch_sh, teasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
