#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"

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
    Debug::Progress progress;

    while (teaReader->ReadEntry()) {
        if (!aaReader->ReadEntry()) {
            Debug(Debug::ERROR) << "AA FASTA has fewer entries than TEA FASTA\n";
            EXIT(EXIT_FAILURE);
        }
        progress.updateProgress();

        const KSeqWrapper::KSeqEntry &teaEntry = teaReader->entry;
        const KSeqWrapper::KSeqEntry &aaEntry = aaReader->entry;

        // Validate same length
        size_t teaLen = teaEntry.sequence.l;
        size_t aaLen = aaEntry.sequence.l;
        if (teaLen != aaLen) {
            Debug(Debug::ERROR) << "Length mismatch at entry " << id
                                << " (TEA header: " << teaEntry.name.s
                                << ", AA header: " << aaEntry.name.s
                                << "): TEA=" << teaLen << " AA=" << aaLen << "\n";
            EXIT(EXIT_FAILURE);
        }

        // Write TEA sequence (as-is, will be interpreted as amino acid characters by submat)
        std::string teaSeq(teaEntry.sequence.s, teaLen);
        teaSeq.push_back('\n');
        teaWriter.writeData(teaSeq.c_str(), teaSeq.length(), id, 0);

        // Write AA sequence
        std::string aaSeq(aaEntry.sequence.s, aaLen);
        aaSeq.push_back('\n');
        aaWriter.writeData(aaSeq.c_str(), aaSeq.length(), id, 0);

        // Write header (use TEA header)
        std::string header(teaEntry.name.s, teaEntry.name.l);
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

    // Check AA FASTA doesn't have more entries
    if (aaReader->ReadEntry()) {
        Debug(Debug::ERROR) << "AA FASTA has more entries than TEA FASTA\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Wrote " << (id - par.identifierOffset) << " TEA+AA sequence pairs\n";

    fclose(lookupFp);
    hdrWriter.close(true);
    aaWriter.close(true);
    teaWriter.close(true);

    delete teaReader;
    delete aaReader;

    return EXIT_SUCCESS;
}
