#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"

#include <string>

// Strip tea_convert entropy suffix ("|H=0.123") from header name
static std::string stripEntropySuffix(const char *name, size_t len) {
    std::string s(name, len);
    size_t pos = s.find("|H=");
    if (pos != std::string::npos) {
        s.erase(pos);
    }
    return s;
}

int createteadb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.printParameters(command.cmd, argc, argv, *command.params);

    // Positional args: <tea_fasta> <aa_fasta> <output_db>
    std::string teaFastaFile = par.filenames[0];
    std::string aaFastaFile = par.filenames[1];
    std::string outDB = par.filenames[2];

    // TEA sequences go into the main DB (used for k-mer prefiltering)
    std::string teaSeqDataFile = outDB;
    std::string teaSeqIndexFile = outDB + ".index";
    std::string hdrDataFile = outDB + "_h";
    std::string hdrIndexFile = outDB + "_h.index";
    std::string lookupFile = outDB + ".lookup";

    // AA sequences go into a companion _aa DB (used for combined alignment scoring)
    std::string aaSeqDataFile = outDB + "_aa";
    std::string aaSeqIndexFile = outDB + "_aa.index";

    DBWriter teaWriter(teaSeqDataFile.c_str(), teaSeqIndexFile.c_str(), 1, 0, Parameters::DBTYPE_AMINO_ACIDS);
    teaWriter.open();
    DBWriter aaWriter(aaSeqDataFile.c_str(), aaSeqIndexFile.c_str(), 1, 0, Parameters::DBTYPE_AMINO_ACIDS);
    aaWriter.open();
    DBWriter hdrWriter(hdrDataFile.c_str(), hdrIndexFile.c_str(), 1, 0, Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();

    FILE *lookupFp = FileUtil::openFileOrDie(lookupFile.c_str(), "w", false);

    KSeqWrapper *teaReader = KSeqFactory(teaFastaFile.c_str());
    KSeqWrapper *aaReader = KSeqFactory(aaFastaFile.c_str());

    unsigned int id = par.identifierOffset;
    unsigned int skipped = 0;
    Debug::Progress progress;

    while (teaReader->ReadEntry()) {
        if (!aaReader->ReadEntry()) {
            Debug(Debug::ERROR) << "AA FASTA has fewer entries than TEA FASTA\n";
            EXIT(EXIT_FAILURE);
        }
        progress.updateProgress();

        const KSeqWrapper::KSeqEntry &teaEntry = teaReader->entry;
        const KSeqWrapper::KSeqEntry &aaEntry = aaReader->entry;

        std::string teaName = stripEntropySuffix(teaEntry.name.s, teaEntry.name.l);
        std::string aaName(aaEntry.name.s, aaEntry.name.l);

        if (teaName != aaName) {
            Debug(Debug::ERROR) << "Header mismatch at entry " << id
                                << ": TEA=" << teaName << " AA=" << aaName
                                << ". TEA and AA FASTA must be in the same order.\n";
            EXIT(EXIT_FAILURE);
        }

        size_t teaLen = teaEntry.sequence.l;
        size_t aaLen = aaEntry.sequence.l;
        if (teaLen != aaLen) {
            Debug(Debug::WARNING) << "Length mismatch for " << teaName
                                  << ": TEA=" << teaLen << " AA=" << aaLen << ", skipping\n";
            skipped++;
            continue;
        }

        // Write TEA sequence
        std::string teaSeq(teaEntry.sequence.s, teaLen);
        teaSeq.push_back('\n');
        teaWriter.writeData(teaSeq.c_str(), teaSeq.length(), id, 0);

        // Write AA sequence
        std::string aaSeq(aaEntry.sequence.s, aaLen);
        aaSeq.push_back('\n');
        aaWriter.writeData(aaSeq.c_str(), aaSeq.length(), id, 0);

        // Write header (use stripped TEA name)
        std::string header = teaName;
        if (teaEntry.comment.l > 0) {
            header.append(" ");
            header.append(teaEntry.comment.s, teaEntry.comment.l);
        }
        header.append("\n");
        hdrWriter.writeData(header.c_str(), header.length(), id, 0);

        // Write lookup entry
        std::string accession = Util::parseFastaHeader(header.c_str());
        fprintf(lookupFp, "%u\t%s\t0\n", id, accession.c_str());

        id++;
    }

    if (aaReader->ReadEntry()) {
        Debug(Debug::ERROR) << "AA FASTA has more entries than TEA FASTA\n";
        EXIT(EXIT_FAILURE);
    }

    unsigned int written = id - par.identifierOffset;
    Debug(Debug::INFO) << "Wrote " << written << " TEA+AA sequence pairs";
    if (skipped > 0) {
        Debug(Debug::INFO) << " (" << skipped << " skipped due to length mismatch)";
    }
    Debug(Debug::INFO) << "\n";

    fclose(lookupFp);
    hdrWriter.close(true);
    aaWriter.close(true);
    teaWriter.close(true);

    delete teaReader;
    delete aaReader;

    return EXIT_SUCCESS;
}
