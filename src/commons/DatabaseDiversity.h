#ifndef STEAM_DATABASE_DIVERSITY_H
#define STEAM_DATABASE_DIVERSITY_H

#include <cstddef>
#include <cstdint>
#include <string>

struct DatabaseDiversityMetadata {
    unsigned int schemaVersion;
    std::string method;
    uint64_t sequences;
    uint64_t residues;
    double effectiveTargets;
    double coordinateStandardDeviation;
    size_t coordinates;
    unsigned int hllPrecision;
};

class DatabaseDiversity {
public:
    static const char *method();
    static const char *suffix();

    // Compute the fixed AA3+TEA5 MinHash/HLL statistic from a paired native
    // STEAM database.  This is a linear, parallel scan and does no clustering.
    static bool compute(const std::string &database, int threads,
                        DatabaseDiversityMetadata &metadata,
                        std::string &error);

    // Metadata writes are atomic.  The reader is strict about method/version
    // so a search cannot silently mix a statistic with incompatible fit
    // coefficients.
    static bool write(const std::string &database,
                      const DatabaseDiversityMetadata &metadata,
                      std::string &error);
    static bool read(const std::string &database,
                     DatabaseDiversityMetadata &metadata,
                     std::string &error);
    static std::string metadataPath(const std::string &database);
    static std::string findMetadataPath(const std::string &database);
};

#endif
