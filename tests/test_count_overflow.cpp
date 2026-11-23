// The count-only prefilter must emit each target once after merging chunks.
#include "CacheFriendlyOperations.h"
#include <cstdio>
#include <map>

template<unsigned Bins> bool check() {
    CacheFriendlyOperations<Bins> counter(64, 32);
    for (unsigned repetition = 0; repetition < 3; ++repetition) {
        CounterResult rows[] = {{7, 5, 9}, {7, 11, 11}, {13, 9, 3}};
        const size_t n = counter.mergeElementsByScore(rows, 3);
        std::map<unsigned, unsigned> observed;
        for (size_t i = 0; i < n; ++i) {
            if (!observed.emplace(unsigned(rows[i].id), unsigned(rows[i].count)).second) {
                std::fprintf(stderr, "Duplicate target %u after count merge (bins=%u)\n", rows[i].id, Bins);
                return false;
            }
        }
        if (n != 2 || observed[7] != 20 || observed[13] != 3) {
            std::fprintf(stderr, "Wrong count merge (bins=%u, n=%zu)\n", Bins, n);
            return false;
        }
        CounterResult saturated[] = {{7, 255, 250}, {7, 17, 10}, {13, 2, 4}};
        const size_t m = counter.mergeElementsByScore(saturated, 3);
        observed.clear();
        for (size_t i = 0; i < m; ++i) {
            if (!observed.emplace(unsigned(saturated[i].id), unsigned(saturated[i].count)).second) return false;
        }
        if (m != 2 || observed[7] != 255 || observed[13] != 4) return false;
    }
    return true;
}
int main() {
    if (!(check<2>() && check<4>() && check<8>() && check<16>() && check<32>() && check<64>())) return 1;
    std::puts("PASS count-overflow merges: unique targets, exact/saturated counts, cleared state");
}
