#include "DatabaseDiversity.h"

#include "DBReader.h"
#include "FileUtil.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <unistd.h>

namespace {

const size_t SKETCHES = 64;
const size_t SKETCH_CAPACITY = 1024;
const uint64_t EMPTY_HASH = std::numeric_limits<uint64_t>::max();
const char *METHOD = "aa3-tea5-doph64-kmv1024-v1";
const char *SUFFIX = ".steam-diversity";

uint64_t splitmix64(uint64_t value) {
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31);
}

char uppercase(char value) {
    return value >= 'a' && value <= 'z' ? value - 32 : value;
}

int symbol(char value) {
    value = uppercase(value);
    switch (value) {
        case 'A': return 0; case 'C': return 1; case 'D': return 2;
        case 'E': return 3; case 'F': return 4; case 'G': return 5;
        case 'H': return 6; case 'I': return 7; case 'K': return 8;
        case 'L': return 9; case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14;
        case 'S': return 15; case 'T': return 16; case 'V': return 17;
        case 'W': return 18; case 'Y': return 19;
        default: return -1;
    }
}

uint64_t contiguousCode(const char *sequence, size_t start, size_t width,
                        bool &valid) {
    uint64_t code = 0;
    valid = true;
    for (size_t position = 0; position < width; ++position) {
        const int value = symbol(sequence[start + position]);
        if (value < 0) {
            valid = false;
            return 0;
        }
        code = code * 20 + static_cast<uint64_t>(value);
    }
    return code;
}

void addFeature(std::array<uint64_t, SKETCHES> &signature, uint64_t feature) {
    const uint64_t hashed = splitmix64(feature);
    const size_t bin = hashed & (SKETCHES - 1);
    signature[bin] = std::min(signature[bin], hashed);
}

std::array<uint64_t, SKETCHES> buildSignature(const char *tea, const char *aa,
                                               size_t length) {
    std::array<uint64_t, SKETCHES> signature;
    signature.fill(EMPTY_HASH);
    if (length >= 3) {
        for (size_t start = 0; start + 3 <= length; ++start) {
            bool valid = false;
            const uint64_t code = contiguousCode(aa, start, 3, valid);
            if (valid) addFeature(signature, code ^ 0xA300000000000000ULL);
        }
    }
    if (length >= 5) {
        for (size_t start = 0; start + 5 <= length; ++start) {
            bool valid = false;
            const uint64_t code = contiguousCode(tea, start, 5, valid);
            if (valid) addFeature(signature, code ^ 0xB500000000000000ULL);
        }
    }

    size_t nonempty = 0;
    for (size_t i = 0; i < SKETCHES; ++i) nonempty += signature[i] != EMPTY_HASH;
    if (nonempty == 0) {
        uint64_t fallback = 0xC700000000000000ULL;
        for (size_t i = 0; i < length; ++i) {
            fallback = splitmix64(
                fallback ^ static_cast<uint8_t>(uppercase(tea[i]))
                ^ (uint64_t(static_cast<uint8_t>(uppercase(aa[i]))) << 8)
            );
        }
        signature.fill(fallback);
        return signature;
    }

    // Densified one-permutation MinHash: deterministically fill empty bins
    // from the next occupied bin without adding a second feature scan.
    for (size_t bin = 0; bin < SKETCHES; ++bin) {
        if (signature[bin] != EMPTY_HASH) continue;
        size_t distance = 1;
        while (signature[(bin + distance) & (SKETCHES - 1)] == EMPTY_HASH) {
            ++distance;
        }
        const uint64_t source = signature[(bin + distance) & (SKETCHES - 1)];
        signature[bin] = splitmix64(
            source ^ (distance * 0xD6E8FEB86659FD93ULL) ^ bin
        );
    }
    return signature;
}

void retainMinimum(std::set<uint64_t> &sketch, uint64_t hash) {
    if (sketch.size() < SKETCH_CAPACITY || hash < *sketch.rbegin()) {
        sketch.insert(hash);
        if (sketch.size() > SKETCH_CAPACITY) sketch.erase(std::prev(sketch.end()));
    }
}

void addSignature(std::set<uint64_t> *sketches,
                  const std::array<uint64_t, SKETCHES> &signature) {
    for (size_t coordinate = 0; coordinate < SKETCHES; ++coordinate) {
        const uint64_t hashed = splitmix64(
            signature[coordinate]
            ^ (coordinate * 0x94D049BB133111EBULL)
        );
        retainMinimum(sketches[coordinate], hashed);
    }
}

template<typename T>
bool readField(std::istream &input, const char *name, T &value) {
    std::string line;
    const std::string prefix = std::string(name) + '\t';
    if (!std::getline(input, line) || line.compare(0, prefix.size(), prefix) != 0)
        return false;
    const std::string text = line.substr(prefix.size());
    if (text.empty() || text[0] == '-' || text[0] == ' ' || text[0] == '\t') return false;
    std::istringstream field(text);
    std::string extra;
    return (field >> value) && !(field >> extra);
}

bool validMetadata(const DatabaseDiversityMetadata &metadata) {
    return metadata.schemaVersion == 2 && metadata.method == METHOD
        && metadata.sequences > 0 && std::isfinite(metadata.effectiveTargets)
        && metadata.effectiveTargets > 0.0
        && std::isfinite(metadata.coordinateStandardDeviation)
        && metadata.coordinateStandardDeviation >= 0.0
        && metadata.coordinates == SKETCHES && metadata.sketchCapacity == SKETCH_CAPACITY;
}

bool hasSuffix(const std::string &value, const std::string &suffix) {
    return value.size() >= suffix.size()
        && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

} // namespace

const char *DatabaseDiversity::method() { return METHOD; }
const char *DatabaseDiversity::suffix() { return SUFFIX; }

std::string DatabaseDiversity::metadataPath(const std::string &database) {
    return database + SUFFIX;
}

std::string DatabaseDiversity::findMetadataPath(const std::string &database) {
    const std::string direct = metadataPath(database);
    if (FileUtil::fileExists(direct.c_str())) return direct;
    for (const char *extension : {".idx", ".linidx"}) {
        const std::string suffixValue(extension);
        if (hasSuffix(database, suffixValue)) {
            const std::string base = metadataPath(
                database.substr(0, database.size() - suffixValue.size())
            );
            if (FileUtil::fileExists(base.c_str())) return base;
        }
    }
    return direct;
}

bool DatabaseDiversity::compute(
        const std::string &database, int threads,
        DatabaseDiversityMetadata &metadata, std::string &error) {
    if (threads < 1) {
        error = "Diversity computation needs at least one thread";
        return false;
    }
    const std::string aaDatabase = database + "_aa";
    DBReader<unsigned int> tea(
        database.c_str(), (database + ".index").c_str(), threads,
        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
    );
    DBReader<unsigned int> aa(
        aaDatabase.c_str(), (aaDatabase + ".index").c_str(), threads,
        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
    );
    tea.open(DBReader<unsigned int>::NOSORT);
    aa.open(DBReader<unsigned int>::NOSORT);
    if (tea.getSize() == 0 || tea.getSize() != aa.getSize()) {
        tea.close();
        aa.close();
        error = "Diversity requires a nonempty paired TEA/AA database";
        return false;
    }
    for (size_t id = 0; id < tea.getSize(); ++id) {
        if (tea.getDbKey(id) != aa.getDbKey(id)
                || tea.getSeqLen(id) != aa.getSeqLen(id)) {
            tea.close();
            aa.close();
            error = "TEA/AA database identities or lengths differ";
            return false;
        }
    }

    int workerCount = threads;
#ifndef _OPENMP
    workerCount = 1;
#endif
    std::vector<std::set<uint64_t>> threadSketches(workerCount * SKETCHES);
    uint64_t residues = 0;

#pragma omp parallel num_threads(workerCount) reduction(+:residues)
    {
#ifdef _OPENMP
        const int thread = omp_get_thread_num();
#else
        const int thread = 0;
#endif
        std::set<uint64_t> *local = threadSketches.data() + thread * SKETCHES;
#pragma omp for schedule(dynamic, 256)
        for (size_t id = 0; id < tea.getSize(); ++id) {
            const size_t length = tea.getSeqLen(id);
            const char *teaSequence = tea.getData(id, thread);
            const char *aaSequence = aa.getData(id, thread);
            addSignature(local, buildSignature(teaSequence, aaSequence, length));
            residues += length;
        }
    }

    std::array<double, SKETCHES> estimates;
#pragma omp parallel for num_threads(workerCount) schedule(static)
    for (size_t coordinate = 0; coordinate < SKETCHES; ++coordinate) {
        std::set<uint64_t> merged;
        for (int thread = 0; thread < workerCount; ++thread) {
            for (uint64_t hash : threadSketches[thread * SKETCHES + coordinate])
                retainMinimum(merged, hash);
        }
        // Exact while unsaturated; otherwise the unbiased bottom-k estimator.
        estimates[coordinate] = merged.size() < SKETCH_CAPACITY ? double(merged.size())
            : double((SKETCH_CAPACITY - 1) * 18446744073709551616.0L
                     / (static_cast<long double>(*merged.rbegin()) + 1));
    }

    double mean = 0.0;
    for (size_t coordinate = 0; coordinate < SKETCHES; ++coordinate) {
        mean += estimates[coordinate] / SKETCHES;
    }
    double variance = 0.0;
    for (size_t coordinate = 0; coordinate < SKETCHES; ++coordinate) {
        const double difference = estimates[coordinate] - mean;
        variance += difference * difference / SKETCHES;
    }
    const uint64_t sequences = static_cast<uint64_t>(tea.getSize());
    tea.close();
    aa.close();
    metadata = {
        2, METHOD, sequences, residues, mean,
        std::sqrt(variance), SKETCHES, SKETCH_CAPACITY
    };
    error.clear();
    return true;
}

bool DatabaseDiversity::write(
        const std::string &database,
        const DatabaseDiversityMetadata &metadata, std::string &error) {
    if (!validMetadata(metadata)) {
        error = "Refusing to write invalid diversity metadata";
        return false;
    }
    const std::string destination = metadataPath(database);
    const std::string temporary = destination + ".tmp." + std::to_string(getpid());
    {
        std::ofstream output(temporary.c_str());
        if (!output) {
            error = "Cannot write diversity metadata: " + temporary;
            return false;
        }
        output << std::setprecision(17)
               << "steam_database_diversity\t" << metadata.schemaVersion << "\n"
               << "method\t" << metadata.method << "\n"
               << "sequences\t" << metadata.sequences << "\n"
               << "residues\t" << metadata.residues << "\n"
               << "effective_targets\t" << metadata.effectiveTargets << "\n"
               << "coordinate_sd\t" << metadata.coordinateStandardDeviation << "\n"
               << "coordinates\t" << metadata.coordinates << "\n"
               << "sketch_capacity\t" << metadata.sketchCapacity << "\n";
        output.close();
        if (!output) {
            std::remove(temporary.c_str());
            error = "Failed writing diversity metadata: " + temporary;
            return false;
        }
    }
    if (std::rename(temporary.c_str(), destination.c_str()) != 0) {
        const std::string message = std::strerror(errno);
        std::remove(temporary.c_str());
        error = "Cannot publish diversity metadata: " + message;
        return false;
    }
    error.clear();
    return true;
}

bool DatabaseDiversity::read(const std::string &database,
                             DatabaseDiversityMetadata &metadata,
                             std::string &error) {
    const std::string path = findMetadataPath(database);
    std::ifstream input(path.c_str());
    if (!input) {
        error = "Missing " + path + "; run `steam computediversity " + database
            + "` or rebuild the database with `steam createdb`";
        return false;
    }
    if (!readField(input, "steam_database_diversity", metadata.schemaVersion)
            || !readField(input, "method", metadata.method)
            || !readField(input, "sequences", metadata.sequences)
            || !readField(input, "residues", metadata.residues)
            || !readField(input, "effective_targets", metadata.effectiveTargets)
            || !readField(input, "coordinate_sd", metadata.coordinateStandardDeviation)
            || !readField(input, "coordinates", metadata.coordinates)
            || !readField(input, "sketch_capacity", metadata.sketchCapacity)
            || input.peek() != std::char_traits<char>::eof() || !validMetadata(metadata)) {
        error = "Invalid or incompatible diversity metadata: " + path
            + "; run `steam computediversity " + database + "`";
        return false;
    }
    error.clear();
    return true;
}
