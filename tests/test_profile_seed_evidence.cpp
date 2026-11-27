// Standalone: c++ -std=c++14 -O2 tests/test_profile_seed_evidence.cpp -o test-profile-seeds
#include "../src/commons/AlignmentSeedEvidence.h"
#include <cassert>
#include <iostream>

static std::string profile(const std::string &sequence) {
    const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
    std::string rows(sequence.size() * 20, char(-100));
    for (size_t i = 0; i < sequence.size(); ++i)
        rows[i * 20 + alphabet.find(sequence[i])] = 10;
    return rows;
}

int main() {
    const char *alphabet = "ACDEFGHIKLMNPQRSTVWY";
    std::string model = "STEAM_SEED_MARKOV_V1\nACDEFGHIKLMNPQRSTVWY\n100\n";
    for (int i = 0; i < 21; ++i) {
        model += "0.7";
        for (int a = 1; a < 20; ++a) model += " 0.015789473684210526";
        model += '\n';
    }
    for (std::string pattern : {"1101101", "1110111"}) {
        AlignmentSeedEvidence engine;
        engine.loadText(model, pattern);
        const int k = std::count(pattern.begin(), pattern.end(), '1');
        for (std::string sequence : {"AAA", "AAAAAAAAA", "ACDEFGHIKLMNPQRSTVWY"}) {
            auto rows = profile(sequence);
            auto q = engine.prepareProfile(rows.data(), sequence.size(), 20, alphabet, 10*k);
            auto exact = engine.prepare(sequence.data(), sequence.size());
            std::vector<std::pair<size_t, std::string>> words;
            for (size_t i = 0; i + 7 <= sequence.size(); ++i)
                words.emplace_back(i, sequence.substr(i, 7));
            auto observed = engine.prepareObserved(sequence.size(), words);
            assert(q.seedPositions == (sequence.size() < 7 ? 0 : sequence.size()-6));
            for (size_t start = 0; start < sequence.size(); ++start) {
                auto target = sequence.substr(start);
                auto trace = std::string(target.size(), 'M');
                double actual = engine.evidence(q, target.data(), target.size(), start, 0, trace);
                // Point profiles share the sequence numerator, but normalize
                // by positions rather than distinct words (different for repeats).
                double expected = q.seedPositions ? engine.evidence(exact, target.data(), target.size(), start, 0, trace)
                    * exact.weights.size() / q.seedPositions : 0.;
                assert(std::fabs(actual - expected) < 1e-12);
                double mapped = engine.evidence(observed, target.data(), target.size(), start, 0, trace);
                assert(std::fabs(mapped - expected) < 1e-12);
                assert(mapped <= AlignmentSeedEvidence::maxEvidence(observed));
                assert(actual <= AlignmentSeedEvidence::maxEvidence(q));
            }
        }
        // Two accepted non-consensus letters per position, with repeated windows.
        auto rows = profile("AAAAAAAAA");
        for (size_t i = 0; i < 9; ++i) rows[i * 20 + 1] = 10;
        auto q = engine.prepareProfile(rows.data(), 9, 20, alphabet, 10*k);
        assert(q.seedPositions == 3);
        const std::string target(9, 'C');
        auto exact = engine.prepare(target.data(), target.size());
        double expected = engine.evidence(exact, target.data(), 9, 0, 0, "MMMMMMMMM") / q.seedPositions;
        assert(engine.evidence(q, target.data(), 9, 0, 0, "MMMMMMMMM") == expected);
        assert(engine.evidence(q, target.data(), 9, 0, 0, "MMMMMMM") == expected);
        assert(engine.evidence(q, target.data(), 9, 0, 0, "MMMIMMM") == 0.);
        assert(engine.evidence(q, target.data(), 9, 0, 0, "MMMDMMM") == 0.);
        assert(engine.evidence(q, "XCCCCCC", 7, 0, 0, "MMMMMMM") == 0.);
        // Brute-force a small alphabet to check threshold acceptance, including
        // a threshold above every possible word score.
        rows = profile("AAAAAAA");
        for (size_t i = 0; i < 7; ++i) rows[i * 20 + 1] = 6;
        for (int threshold = 6*k; threshold <= 10*k+1; ++threshold) {
            q = engine.prepareProfile(rows.data(), 7, 20, alphabet, threshold);
            for (int word = 0; word < (1 << k); ++word) {
                int score = 10*k - 4*__builtin_popcount(unsigned(word));
                std::string sequence(7, 'A');
                int bit = 0;
                for (size_t pos = 0; pos < pattern.size(); ++pos)
                    if (pattern[pos] == '1') sequence[pos] = (word & (1 << bit++)) ? 'C' : 'A';
                auto point = engine.prepare(sequence.data(), sequence.size());
                double expected = score >= threshold
                    ? engine.evidence(point, sequence.data(), 7, 0, 0, "MMMMMMM") : 0.;
                assert(engine.evidence(q, sequence.data(), 7, 0, 0, "MMMMMMM") == expected);
            }
        }
        // Fully diffuse rows have the same position denominator. Repeated
        // supported words still count only once, and the early-reject bound holds.
        rows.assign(9 * 20, 0);
        q = engine.prepareProfile(rows.data(), 9, 20, alphabet, 0);
        assert(q.seedPositions == 3);
        assert(engine.evidence(q, target.data(), 9, 0, 0, "MMMMMMMMM") == expected);
        for (const std::string sequence : {target, std::string("ACDEFGHIK")})
            assert(engine.evidence(q, sequence.data(), 9, 0, 0, "MMMMMMMMM") <= AlignmentSeedEvidence::maxEvidence(q));
        // Membership is position-specific, not inferred from permissive rows.
        // Duplicates at one or several starts contribute only once per alignment.
        auto observed = engine.prepareObserved(9, {{0, "CCCCCCC"}, {2, "CCCCCCC"}, {0, "CCCCCCC"}});
        assert(observed.seedPositions == 3 && observed.weights.size() == 1);
        const double rarity = engine.evidence(exact, target.data(), 9, 0, 0, "MMMMMMMMM");
        assert(engine.evidence(observed, target.data(), 9, 0, 0, "MMMMMMMMM") == rarity / 3.);
        assert(engine.evidence(observed, target.data(), 9, 1, 0, "MMMMMMM") == 0.);
        assert(engine.evidence(observed, target.data(), 9, 2, 0, "MMMMMMM") == rarity / 3.);
        for (const std::string trace : {"MMMIMMM", "MMMDMMM"})
            assert(engine.evidence(observed, target.data(), 9, 0, 0, trace) == 0.);
        assert(engine.evidence(observed, "AAAAAAA", 7, 0, 0, "MMMMMMM") == 0.);
        std::string spaced(7, 'C');
        for (size_t i = 0; i < pattern.size(); ++i)
            if (pattern[i] == '0') spaced[i] = 'A';
        assert(engine.evidence(observed, spaced.data(), 7, 0, 0, "MMMMMMM") == rarity / 3.);
        assert(AlignmentSeedEvidence::maxEvidence(observed) >= rarity);
        for (size_t length : {size_t(3), size_t(9)}) {
            auto empty = engine.prepareObserved(length, {});
            assert(AlignmentSeedEvidence::maxEvidence(empty) == 0.);
            assert(engine.evidence(empty, target.data(), 9, 0, 0, std::string(length, 'M')) == 0.);
        }
        assert(AlignmentSeedEvidence::score(100, rarity / 3.) == int(std::floor(100. * (1. + rarity / 3.) + .5)));
        std::cout << "PASS profile and observed seed evidence " << pattern << '\n';
    }
}
