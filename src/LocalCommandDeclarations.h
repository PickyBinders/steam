#ifndef STEAM_LOCAL_COMMAND_DECLARATIONS_H
#define STEAM_LOCAL_COMMAND_DECLARATIONS_H

#include "Command.h"

extern int createteadb(int argc, const char **argv, const Command& command);
extern int teaalign(int argc, const char **argv, const Command& command);
extern int steamconvertalis(int argc, const char **argv, const Command& command);
extern int teasearch(int argc, const char **argv, const Command& command);
extern int easyteasearch(int argc, const char **argv, const Command& command);
extern int teacluster(int argc, const char **argv, const Command& command);
extern int easyteacluster(int argc, const char **argv, const Command& command);

#endif // STEAM_LOCAL_COMMAND_DECLARATIONS_H
