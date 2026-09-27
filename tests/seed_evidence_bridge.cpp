// Standalone diagnostic bridge: model + pattern arguments, whitespace rows
// query target qstart0 tstart0 uncompressed_backtrace raw_score on stdin.
#include "../src/commons/AlignmentSeedEvidence.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char **argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: seed_evidence_bridge MODEL PATTERN");
        AlignmentSeedEvidence model;
        model.load(argv[1], argv[2]);
        std::string q, t, bt;
        size_t qs, ts;
        int raw;
        std::cout << std::setprecision(17);
        while (std::cin >> q >> t >> qs >> ts >> bt >> raw) {
            if (q == "-") q.clear();
            if (t == "-") t.clear();
            if (bt == "-") bt.clear();
            const auto prepared = model.prepare(q.data(), q.size());
            double e = model.evidence(prepared, t.data(), t.size(), qs, ts, bt);
            std::cout << e << '\t' << AlignmentSeedEvidence::score(raw, e)
                      << '\t' << prepared.weights.size() << '\n';
        }
        if (!std::cin.eof()) throw std::runtime_error("Malformed bridge input");
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
