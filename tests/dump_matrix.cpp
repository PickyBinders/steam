// Inspect the native loaded scores for comparison with TEA's integer matrix.
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include <cstdio>
#include <cstdlib>

const char* binary_name = "steam_dump_matrix";
DEFAULT_PARAMETER_SINGLETON_INIT

int main(int argc, const char** argv) {
    if (argc != 3) return 2;
    SubstitutionMatrix matrix(argv[1], std::atof(argv[2]), 0);
    for (int i = 0; i < matrix.alphabetSize; ++i)
        for (int j = 0; j < matrix.alphabetSize; ++j)
            std::printf("%c\t%c\t%d\n", matrix.num2aa[i], matrix.num2aa[j], matrix.subMatrix[i][j]);
    return 0;
}
