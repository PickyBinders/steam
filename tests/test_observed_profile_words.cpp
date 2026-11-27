// c++ -std=c++14 -Ilib/mmseqs/src/prefiltering tests/test_observed_profile_words.cpp -o test-words
#include "ObservedProfileWords.h"
#include <cassert>
#include <sstream>

int main(int argc, char **argv) {
    ObservedProfileWords words;
    if (argc > 1) {
        std::istringstream invalid(argv[1]);
        words.read(invalid);
        words.get(0, 9);
        return 0;
    }
    std::istringstream input("0 0 ACDAFGH\n0 2 DEFAHIK\n0 0 ACDAFGH\n7 0 CCCACCC\n");
    words.read(input);
    const auto &first = words.get(0, 9);
    assert(first.size() == 3 && first[1].first == 2 && first[1].second == "DEFAHIK");
    assert(words.get(7, 7).size() == 1);
    assert(words.get(1, 3).empty());
    std::istringstream empty("");
    words.read(empty);
    assert(words.get(0, 9).empty());
}
