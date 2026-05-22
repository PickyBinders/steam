#ifndef STEAM_LOCAL_COMMAND_DECLARATIONS_H
#define STEAM_LOCAL_COMMAND_DECLARATIONS_H

#include "Command.h"

// steam-native commands
extern int createteadb(int argc, const char **argv, const Command& command);
extern int createteasubdb(int argc, const char **argv, const Command& command);
extern int teaalign(int argc, const char **argv, const Command& command);
extern int tearescorediagonal(int argc, const char **argv, const Command& command);
extern int steamconvertalis(int argc, const char **argv, const Command& command);
extern int teasearch(int argc, const char **argv, const Command& command);
extern int easyteasearch(int argc, const char **argv, const Command& command);
extern int prefilter(int argc, const char **argv, const Command& command);

// base mmseqs commands we wrap to set env vars + delegate
extern int clusteringworkflow(int argc, const char **argv, const Command& command);
extern int linclust(int argc, const char **argv, const Command& command);
extern int easycluster(int argc, const char **argv, const Command& command);
extern int easylinclust(int argc, const char **argv, const Command& command);

// steam wrappers (defined in src/workflow/TeaClusterWrappers.cpp)
extern int teacluster(int argc, const char **argv, const Command& command);
extern int teaLinclust(int argc, const char **argv, const Command& command);
extern int easyteacluster(int argc, const char **argv, const Command& command);
extern int easytealinclust(int argc, const char **argv, const Command& command);

#endif // STEAM_LOCAL_COMMAND_DECLARATIONS_H
