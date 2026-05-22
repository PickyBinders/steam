#include "Command.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"
#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"

LocalParameters& localPar = LocalParameters::getLocalInstance();

static const char* STEAM_AUTHOR = "Janani Durairaj <janani.durairaj@unibas.ch>";

// Required by MMseqs2 framework
void updateValidation() {}
void (*validatorUpdate)(void) = updateValidation;

std::vector<Command> steamCommands = {
        // Easy workflows
        {"easy-search",  easyteasearch,  &localPar.easyteasearchworkflow,  COMMAND_EASY,
                "Sensitive homology search with combined TEA+AA scoring",
                "# Search FASTA pairs against a database\n"
                "steam easy-search query_tea.fa query_aa.fa target_tea.fa target_aa.fa result.m8 tmp\n\n"
                "# Search FASTA pair against pre-built DB\n"
                "steam easy-search query_tea.fa query_aa.fa targetDB result.m8 tmp\n",
                STEAM_AUTHOR,
                "<i:queryTeaFasta> <i:queryAAFasta> <i:targetTeaFasta> <i:targetAAFasta>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_MMSEQS2, {{"inputFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                   {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-cluster",         easyteacluster,    &localPar.easyclusterworkflow,    COMMAND_EASY,
                "Cluster paired TEA/AA FASTAs with cascaded TEA-aware scoring",
                "steam easy-cluster tea.fasta aa.fasta result tmp\n",
                STEAM_AUTHOR,
                "<i:teaFasta> <i:aaFasta> <o:clusterPrefix> <tmpDir>",
                CITATION_MMSEQS2, {{"inputFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                   {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-linclust",        easytealinclust,   &localPar.easylinclustworkflow,   COMMAND_EASY,
                "Linear-time clustering of paired TEA/AA FASTAs with TEA-aware scoring",
                "steam easy-linclust tea.fasta aa.fasta result tmp\n",
                STEAM_AUTHOR,
                "<i:teaFasta> <i:aaFasta> <o:clusterPrefix> <tmpDir>",
                CITATION_MMSEQS2, {{"inputFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                   {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        // Database creation
        {"createdb",     createteadb,  &localPar.createteadb,  COMMAND_DATABASE_CREATION,
                "Convert paired TEA/AA FASTA files to a sequence DB",
                "# Create a database from paired FASTA files\n"
                "steam createdb tea.fasta aa.fasta steamDB\n",
                STEAM_AUTHOR,
                "<i:teaFastaFile> <i:aaFastaFile> <o:steamDB>",
                CITATION_MMSEQS2, {{"teaFastaFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"aaFastaFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"steamDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        // createsubdb shadow: also subsets the _aa companion DB
        {"createsubdb",  createteasubdb, &localPar.createsubdb, COMMAND_SET,
                "createsubdb that also subsets the steam _aa companion DB",
                "steam createsubdb subsetList inDB outDB\n",
                STEAM_AUTHOR,
                "<i:subsetFile|DB> <i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"subsetFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDbAndFlat },
                                   {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                   {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        // Core search workflow (DB input/output)
        {"search",       teasearch,    &localPar.teasearchworkflow,     COMMAND_MAIN,
                "Sensitive homology search with combined TEA+AA scoring",
                "steam search queryDB targetDB result tmp\n",
                STEAM_AUTHOR,
                "<i:queryDB> <i:targetDB> <o:resultDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"cluster",      teacluster,   &localPar.clusterworkflow,       COMMAND_MAIN,
                "Cascaded clustering with TEA-aware scoring",
                "steam cluster sequenceDB clusterDB tmp\n",
                STEAM_AUTHOR,
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"linclust",     teaLinclust,  &localPar.linclustworkflow,      COMMAND_MAIN,
                "Linear-time clustering with TEA-aware scoring",
                "steam linclust sequenceDB clusterDB tmp\n",
                STEAM_AUTHOR,
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        // Prefilter (MMseqs2 k-mer prefilter, on TEA sequences with steam defaults)
        {"prefilter",    prefilter,    &localPar.teaprefilter,            COMMAND_PREFILTER,
                "Fast k-mer prefilter on TEA sequences",
                "steam prefilter queryDB targetDB prefilterDB\n",
                STEAM_AUTHOR,
                "<i:queryDB> <i:targetDB> <o:prefilterDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefilterDb }}},
        // rescorediagonal shadow: routes mmseqs cluster's hamming/ungapped steps through TEA scoring
        {"rescorediagonal", tearescorediagonal, &localPar.tearescorediagonal, COMMAND_ALIGNMENT,
                "TEA-aware rescorediagonal (shadow over base rescorediagonal)",
                "steam rescorediagonal queryDB targetDB prefiltResult rescoreResult\n",
                STEAM_AUTHOR,
                "<i:queryDB> <i:targetDB> <i:prefiltResult> <o:rescoreResult>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"prefiltResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"rescoreResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        // Alignment
        {"align",        teaalign,     &localPar.teaalign,              COMMAND_ALIGNMENT,
                "Align with combined TEA+AA dual scoring",
                "steam align queryDB targetDB prefiltResult alnResult\n",
                STEAM_AUTHOR,
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
                STEAM_AUTHOR,
                "<i:queryDB> <i:targetDB> <i:alnResult> <o:outputFile>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                   {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                   {"outputFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
};

// Empty external thresholds and downloads (required by MMseqs2 framework)
std::vector<KmerThreshold> externalThreshold = {};
std::vector<DatabaseDownload> externalDownloads = {};
