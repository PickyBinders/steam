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

namespace {

constexpr const char *STEAM_K5_PATTERN = "1101101";
constexpr const char *STEAM_K6_PATTERN = "1110111";

void resolveSteamKmerPattern(LocalParameters &par) {
    if (par.PARAM_SPACED_KMER_PATTERN.wasSet && par.kmerSize == 0) {
        par.kmerSize = 0;
        for (char position : par.spacedKmerPattern) {
            par.kmerSize += position == '1';
        }
    }
    // Read only target index metadata, never target sequence data. Scale the
    // MMseqs k6 limit by one 20-letter seed position for STEAM's k5/k6 range.
    if (par.kmerSize == 0) {
        IndexReader target(par.filenames[1], 1, IndexReader::SEQUENCES,
                           IndexReader::PRELOAD_NO, DBReader<unsigned int>::USE_INDEX);
        if (target.index != NULL) {
            par.kmerSize = PrefilteringIndexReader::getMetadata(target.index).kmerSize;
            if (!par.PARAM_SPACED_KMER_PATTERN.wasSet) {
                par.spacedKmerPattern = PrefilteringIndexReader::getSpacedPattern(target.index);
                par.spacedKmer = PrefilteringIndexReader::getMetadata(target.index).spacedKmer;
            }
        } else {
            const size_t k5Limit = IndexTable::getUpperBoundAACountForKmerSize(6) / 20;
            par.kmerSize = target.sequenceReader->getAminoAcidDBSize() < k5Limit ? 5 : 6;
        }
    }
    if (par.spacedKmer && !par.PARAM_SPACED_KMER_PATTERN.wasSet
            && par.spacedKmerPattern.empty()) {
        if (par.kmerSize == 5) {
            par.spacedKmerPattern = STEAM_K5_PATTERN;
        } else if (par.kmerSize == 6) {
            par.spacedKmerPattern = STEAM_K6_PATTERN;
        }
    }
    if (!par.spacedKmerPattern.empty()) {
        int patternWeight = 0;
        for (char position : par.spacedKmerPattern) {
            patternWeight += position == '1';
        }
        if (patternWeight != par.kmerSize) {
            Debug(Debug::ERROR)
                << "STEAM spaced k-mer pattern weight " << patternWeight
                << " does not match k=" << par.kmerSize << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO)
        << "STEAM k-mer selection: k=" << par.kmerSize
        << ", pattern=" << (par.spacedKmerPattern.empty()
                              ? "<contiguous>" : par.spacedKmerPattern)
        << "\n";
}

} // namespace

int teasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    teaSearchDefault(par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.teaMatrixFile.empty()) {
        Debug(Debug::ERROR) << "--matcha is required for steam search\n";
        EXIT(EXIT_FAILURE);
    }
    par.validateCalibration();
    resolveSteamKmerPattern(par);
    if (!par.PARAM_SEED_PATTERN.wasSet) par.seedPattern = par.spacedKmerPattern;
    if (par.seedCorrection && par.seedPattern != par.spacedKmerPattern) {
        Debug(Debug::ERROR) << "Alignment and prefilter seed patterns must agree\n";
        EXIT(EXIT_FAILURE);
    }
    // Exhaustive alignment bypasses candidate selection, not final score terms.
    const bool combinedUngapped = par.ungappedTeaAa && !par.exhaustiveSearch;

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
    const float alignmentCompBiasCorrectionScale = par.compBiasCorrectionScale;
    par.scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>(par.teaMatrixFile, par.teaMatrixFile));
    par.compBiasCorrectionScale = 0.15;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.teaprefilter).c_str());
    cmd.addVariable("STEAM_SPLIT_INVARIANT", "1");
    cmd.addVariable("STEAM_UNGAPPED_AA", combinedUngapped ? "1" : NULL);
    cmd.addVariable("STEAM_UNGAPPED_AA_SUBMAT",
                    combinedUngapped ? origScoringMatrixFile.values.aminoacid().c_str() : NULL);
    const std::string aaWeight = SSTR(par.teaWeight);
    cmd.addVariable("STEAM_UNGAPPED_AA_WEIGHT", combinedUngapped ? aaWeight.c_str() : NULL);

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
