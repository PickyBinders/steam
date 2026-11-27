#ifndef STEAM_ALIGNMENT_SEED_EVIDENCE_H
#define STEAM_ALIGNMENT_SEED_EVIDENCE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Alignment-supported S*(1+e), applied once after traceback. Significance
// is calibrated separately on the corrected score.
class AlignmentSeedEvidence {
public:
    enum { SPAN = 7, LETTERS = 20 };
    struct Query {
        size_t length = 0;
        std::vector<int> codes;
        std::unordered_map<int, double> weights;
        // Profile rows use the Markov alphabet order and native eighth-bit scores.
        bool profile = false;
        int threshold = 0;
        std::vector<std::array<int, LETTERS>> rows;
        bool observed = false;
        std::vector<std::vector<int>> observedWords;
        size_t seedPositions = 0;
        double maximum = 0.;
    };

    void load(const std::string &path, const std::string &pattern) {
        std::ifstream in(path.c_str());
        read(in, pattern);
    }

    void loadText(const std::string &text, const std::string &pattern) {
        std::istringstream in(text);
        read(in, pattern);
    }

    void read(std::istream &in, const std::string &pattern) {
        if (pattern != "1101101" && pattern != "1110111")
            fail("Seed evidence requires explicit pattern 1101101 (k5) or 1110111 (k6)");
        offsets.clear();
        for (size_t i = 0; i < pattern.size(); ++i)
            if (pattern[i] == '1') offsets.push_back(i);
        std::string magic, alphabet, trailing;
        double windows;
        if (!(in >> magic >> alphabet >> windows) || magic != "STEAM_SEED_MARKOV_V1"
                || alphabet.size() != LETTERS || !positive(windows))
            fail("Invalid seed Markov model header");
        std::string sorted = alphabet;
        std::sort(sorted.begin(), sorted.end());
        if (sorted != "ACDEFGHIKLMNPQRSTVWY") fail("Invalid seed alphabet");
        code.fill(-1);
        for (int i = 0; i < LETTERS; ++i) {
            const unsigned char c = alphabet[i];
            code[c] = code[c + ('a' - 'A')] = i;
        }
        double sum = 0;
        for (int i = 0; i < LETTERS; ++i) {
            double p;
            if (!(in >> p) || !positive(p)) fail("Invalid seed probabilities");
            sum += p; logP[i] = std::log(p);
        }
        if (std::fabs(sum - 1.) > 1e-8) fail("Seed probabilities must sum to one");
        double transition[LETTERS][LETTERS];
        for (int i = 0; i < LETTERS; ++i) {
            sum = 0;
            for (int j = 0; j < LETTERS; ++j) {
                double &p = transition[i][j];
                if (!(in >> p) || !positive(p)) fail("Invalid seed transitions");
                sum += p;
                logT[0][i][j] = std::log(p);
            }
            if (std::fabs(sum - 1.) > 1e-8) fail("Seed transition row must sum to one");
        }
        if (in >> trailing) fail("Unexpected trailing seed model fields");
        for (int i = 0; i < LETTERS; ++i)
            for (int j = 0; j < LETTERS; ++j) {
                double p = 0;
                for (int k = 0; k < LETTERS; ++k) p += transition[i][k] * transition[k][j];
                logT[1][i][j] = std::log(p);
            }
        logWindows = std::log(windows);
    }

    Query prepare(const char *sequence, size_t length) const {
        Query q;
        q.length = length;
        if (length < SPAN) return q;
        q.codes.reserve(length - SPAN + 1);
        for (size_t i = 0; i + SPAN <= length; ++i) {
            const int word = pack(sequence + i);
            q.codes.push_back(word);
            if (word >= 0 && q.weights.find(word) == q.weights.end()) {
                int previous = code[static_cast<unsigned char>(sequence[i + offsets[0]])];
                double logProbability = logP[previous];
                for (size_t j = 1; j < offsets.size(); ++j) {
                    int next = code[static_cast<unsigned char>(sequence[i + offsets[j]])];
                    logProbability += logT[offsets[j] - offsets[j-1] - 1][previous][next];
                    previous = next;
                }
                q.weights.emplace(word, std::max(0., -logProbability - logWindows));
            }
        }
        return q;
    }

    // A traceback can realize only one target word per profile start position.
    // Normalize by these opportunities, not every alternative accepted word.
    Query prepareProfile(const char *rows, size_t length, size_t stride,
                         const char *alphabet, int threshold) const {
        Query q;
        q.length = length;
        q.profile = true;
        q.threshold = std::max(0, threshold);
        q.rows.resize(length);
        std::array<bool, LETTERS> seen = {};
        if (stride < LETTERS) fail("Invalid seed profile stride");
        for (size_t a = 0; a < LETTERS; ++a) {
            int letter = code[static_cast<unsigned char>(alphabet[a])];
            if (letter < 0 || seen[letter]) fail("Invalid seed profile alphabet");
            seen[letter] = true;
            for (size_t i = 0; i < length; ++i)
                q.rows[i][letter] = static_cast<int8_t>(rows[i * stride + a]);
        }
        if (length < SPAN) return q;
        q.seedPositions = length - SPAN + 1;
        // A conservative bound: every supported word's rarity is at most this.
        double logMinimum = *std::min_element(logP, logP + LETTERS);
        for (size_t j = 1; j < offsets.size(); ++j) {
            const auto &transition = logT[offsets[j] - offsets[j-1] - 1];
            double lowest = 0.;
            for (const auto &row : transition)
                lowest = std::min(lowest, *std::min_element(row, row + LETTERS));
            logMinimum += lowest;
        }
        q.maximum = std::max(0., -logMinimum - logWindows) * (1. + 1e-12);
        return q;
    }

    // Exact words mapped from training members onto the final joint profile.
    // Keep the profile-position denominator and deduplicate supported word keys.
    Query prepareObserved(size_t length,
                          const std::vector<std::pair<size_t, std::string>> &words) const {
        Query q;
        q.length = length;
        q.profile = q.observed = true;
        q.seedPositions = length < SPAN ? 0 : length - SPAN + 1;
        q.observedWords.resize(q.seedPositions);
        for (const auto &entry : words) {
            if (entry.first >= q.seedPositions || entry.second.size() != SPAN)
                fail("Observed seed outside profile");
            int word = pack(entry.second.data());
            if (word < 0) fail("Invalid observed seed symbol");
            q.observedWords[entry.first].push_back(word);
            double weight = rarity(word);
            q.weights.emplace(word, weight);
            q.maximum = std::max(q.maximum, weight);
        }
        for (auto &row : q.observedWords) {
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
        }
        q.maximum *= 1. + 1e-12;
        return q;
    }

    double evidence(const Query &q, const char *target, size_t targetLength,
                    size_t qi, size_t ti, const std::string &backtrace) const {
        if (backtrace.empty()) return 0.;
        if (qi > q.length || ti > targetLength)
            fail("Seed-evidence start outside sequence");
        std::vector<int> supported;
        size_t run = 0;
        for (char op : backtrace) {
            if (op != 'M' && op != 'I' && op != 'D')
                fail("Unexpected seed-evidence traceback operation");
            if (op != 'D') ++qi;
            if (op != 'I') ++ti;
            if (qi > q.length || ti > targetLength)
                fail("Seed-evidence traceback outside sequence");
            run = op == 'M' ? run + 1 : 0;
            if (run >= SPAN) {
                const char *seed = target + ti - SPAN;
                if (q.profile) {
                    int word = pack(seed);
                    if (word < 0) continue;
                    if (q.observed) {
                        const auto &words = q.observedWords[qi - SPAN];
                        if (std::binary_search(words.begin(), words.end(), word))
                            supported.push_back(word);
                    } else {
                        int score = 0;
                        for (size_t offset : offsets)
                            score += q.rows[qi - SPAN + offset][code[static_cast<unsigned char>(seed[offset])]];
                        if (score >= q.threshold) supported.push_back(word);
                    }
                } else {
                    int word = q.codes[qi - SPAN];
                    if (word >= 0 && word == pack(seed)) supported.push_back(word);
                }
            }
        }
        std::sort(supported.begin(), supported.end());
        supported.erase(std::unique(supported.begin(), supported.end()), supported.end());
        double total = 0;
        for (int word : supported)
            total += q.profile && !q.observed ? rarity(word) : q.weights.find(word)->second;
        const size_t denominator = q.profile ? q.seedPositions : q.weights.size();
        return denominator ? total / denominator : 0.;
    }

    // A profile has at most seedPositions supported words, each <= maximum.
    // For sequences, match evidence()'s summation order to bound rounding too.
    static double maxEvidence(const Query &q) {
        if (q.profile) return q.seedPositions ? q.maximum : 0.;
        if (q.weights.empty()) return 0.;
        std::vector<int> words;
        words.reserve(q.weights.size());
        for (const auto &word : q.weights) words.push_back(word.first);
        std::sort(words.begin(), words.end());
        double total = 0.;
        for (int word : words) total += q.weights.find(word)->second;
        return total / q.weights.size();
    }

    static int score(int raw, double evidence) {
        const double rounded = std::floor(static_cast<double>(raw) * (1. + evidence) + .5);
        if (raw < 0 || evidence < 0 || !std::isfinite(evidence) || !std::isfinite(rounded)
                || rounded > std::numeric_limits<int>::max())
            fail("Seed score outside signed integer range");
        return static_cast<int>(rounded);
    }

private:
    // STEAM/MMseqs is built with exceptions disabled; fail loudly like its
    // other input guards instead of changing the engine's compiler settings.
    static void fail(const char *message) {
        std::fprintf(stderr, "%s\n", message);
        std::exit(EXIT_FAILURE);
    }
    static bool positive(double p) { return std::isfinite(p) && p > 0.; }
    double rarity(int word) const {
        std::array<int, SPAN> letters;
        for (size_t j = offsets.size(); j-- > 0;) {
            letters[j] = word % LETTERS;
            word /= LETTERS;
        }
        double logProbability = logP[letters[0]];
        for (size_t j = 1; j < offsets.size(); ++j)
            logProbability += logT[offsets[j] - offsets[j-1] - 1][letters[j-1]][letters[j]];
        return std::max(0., -logProbability - logWindows);
    }
    int pack(const char *sequence) const {
        int word = 0;
        for (size_t offset : offsets) {
            int letter = code[static_cast<unsigned char>(sequence[offset])];
            if (letter < 0) return -1;
            word = word * LETTERS + letter;
        }
        return word;
    }
    std::array<int, 256> code;
    std::vector<size_t> offsets;
    double logP[LETTERS], logT[2][LETTERS][LETTERS], logWindows = 0;
};
#endif
