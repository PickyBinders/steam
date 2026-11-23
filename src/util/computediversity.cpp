#include "DatabaseDiversity.h"
#include "Debug.h"
#include "LocalParameters.h"

#include <string>

int computediversity(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    const std::string database = par.filenames[0];
    DatabaseDiversityMetadata metadata;
    std::string error;
    if (!DatabaseDiversity::compute(database, par.threads, metadata, error)
            || !DatabaseDiversity::write(database, metadata, error)) {
        Debug(Debug::ERROR) << "Database diversity failed: " << error << "\n";
        EXIT(EXIT_FAILURE);
    }
    Debug(Debug::INFO)
        << "STEAM database diversity: method=" << metadata.method
        << " sequences=" << metadata.sequences
        << " effective_targets=" << metadata.effectiveTargets
        << " metadata=" << DatabaseDiversity::metadataPath(database) << "\n";
    return EXIT_SUCCESS;
}
