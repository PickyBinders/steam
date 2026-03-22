#include "Command.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"
#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"

LocalParameters& localPar = LocalParameters::getLocalInstance();

// Required by MMseqs2 framework
void updateValidation() {}
void (*validatorUpdate)(void) = updateValidation;

std::vector<Command> steamCommands = {
        // Easy workflows
        {"easy-search",  easyteasearch,  &localPar.easyteasearchworkflow,  COMMAND_EASY,
                "Sensitive homology search with combined TEA+AA scoring",
                "# Search FASTA pairs against a database\n"
                "steam easy-search query_tea.fa query_aa.fa target_tea.fa target_aa.fa result.m8 tmp --sub-mat matcha.out\n\n"
                "# Search FASTA pair against pre-built DB\n"
                "steam easy-search query_tea.fa query_aa.fa targetDB result.m8 tmp --sub-mat matcha.out\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:queryTeaFasta> <i:queryAAFasta> <i:targetTeaFasta> <i:targetAAFasta>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_MMSEQS2, {{"inputFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                   {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-cluster",  easyteacluster,  &localPar.easyteaclusterworkflow,  COMMAND_EASY,
                "Cluster sequences with combined TEA+AA scoring",
                "# Cluster TEA+AA sequences\n"
                "steam easy-cluster tea.fasta aa.fasta result tmp --sub-mat matcha.out\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:teaFasta> <i:aaFasta> <o:clusterPrefix> <tmpDir>",
                CITATION_MMSEQS2, {{"teaFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"aaFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        // Database creation
        {"createdb",     createteadb,  &localPar.createteadb,  COMMAND_DATABASE_CREATION,
                "Convert paired TEA/AA FASTA files to a sequence DB",
                "# Create a database from paired FASTA files\n"
                "steam createdb tea.fasta aa.fasta steamDB\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:teaFastaFile> <i:aaFastaFile> <o:steamDB>",
                CITATION_MMSEQS2, {{"teaFastaFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"aaFastaFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"steamDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        // Core search/cluster workflows (DB input/output)
        {"search",       teasearch,    &localPar.teasearchworkflow,     COMMAND_MAIN,
                "Sensitive homology search with combined TEA+AA scoring",
                "steam search queryDB targetDB result tmp --sub-mat matcha.out\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:queryDB> <i:targetDB> <o:resultDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"cluster",      teacluster,   &localPar.teaclusterworkflow,    COMMAND_MAIN,
                "Cluster sequences with combined TEA+AA scoring",
                "steam cluster steamDB result tmp --sub-mat matcha.out\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:steamDB> <o:clusterDB> <tmpDir>",
                CITATION_MMSEQS2, {{"steamDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        // Ungapped diagonal rescoring
        {"tearescorediagonal", tearescorediagonal, &localPar.tearescorediagonal, COMMAND_ALIGNMENT,
                "Rescore prefilter diagonals with TEA+AA combined ungapped scoring",
                "steam tearescorediagonal queryDB targetDB prefiltResult rescoreResult --tea-mat matcha.out\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:queryDB> <i:targetDB> <i:prefiltResult> <o:rescoreResult>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"prefiltResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"rescoreResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        // Alignment
        {"align",        teaalign,     &localPar.teaalign,              COMMAND_ALIGNMENT,
                "Align with combined TEA+AA dual scoring",
                "steam align queryDB targetDB prefiltResult alnResult --sub-mat matcha.out --tea-weight 1.4\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:queryDB> <i:targetDB> <i:prefiltResult> <o:alnResult>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"prefiltResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"alnResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        // Format conversion
        {"convertalis",  steamconvertalis,  &localPar.convertalignments,  COMMAND_FORMAT_CONVERSION,
                "Convert alignment DB to BLAST-tab, SAM or custom format (with TEA columns)",
                "# Convert to BLAST-tab with TEA identity\n"
                "steam convertalis queryDB targetDB alnResult result.m8 --format-output query,target,fident,tfident,evalue\n\n"
                "# TEA-specific output columns: tfident, tpident, qteaseq, tteaseq, qteaaln, tteaaln\n"
                "# fident/pident = amino acid identity, tfident/tpident = TEA identity\n",
                "Janani Durairaj <janani.durairaj@unibas.ch>",
                "<i:queryDB> <i:targetDB> <i:alnResult> <o:outputFile>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"outputFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
};

// Empty external thresholds and downloads (required by MMseqs2 framework)
std::vector<KmerThreshold> externalThreshold = {};
std::vector<DatabaseDownload> externalDownloads = {};
