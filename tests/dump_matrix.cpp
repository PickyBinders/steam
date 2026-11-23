// Inspect the native loaded scores for comparison with TEA's integer matrix.
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

const char* binary_name = "steam_dump_matrix";
DEFAULT_PARAMETER_SINGLETON_INIT

int main(int argc, const char** argv) {
    if (argc != 3 && argc != 4) return 2;
    SubstitutionMatrix matrix(argv[1], std::atof(argv[2]), 0);
    for (int i = 0; i < matrix.alphabetSize; ++i)
        for (int j = 0; j < matrix.alphabetSize; ++j)
            std::printf("%c\t%c\t%d\n", matrix.num2aa[i], matrix.num2aa[j], matrix.subMatrix[i][j]);
    if (argc == 4) {
        const int length=std::strlen(argv[3]);
        std::vector<unsigned char> sequence(length);
        std::vector<float> bias(length);
        for (int i=0; i<length; ++i) sequence[i]=matrix.aa2num[static_cast<unsigned char>(argv[3][i])];
        SubstitutionMatrix::calcLocalAaBiasCorrection(&matrix,sequence.data(),length,bias.data(),1.0f);
        for (int i=0; i<length; ++i)
            std::printf("BIAS\t%d\t%d\n",i,static_cast<int>(bias[i]<0 ? bias[i]-0.5f : bias[i]+0.5f));
    }
    return 0;
}
