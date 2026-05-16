#include <cstddef>
#include "Command.h"
#include "LocalParameters.h"

const char* binary_name = "steam";
const char* tool_name = "steam";
const char* tool_introduction = "STEAM (Search with TEA against Many)";
const char* main_author = "Janani Durairaj";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
const char* index_version_compatible = "st1";
bool hide_base_commands = true;
bool hide_base_downloads = true;

extern std::vector<Command> baseCommands;
extern std::vector<Command> steamCommands;
void init() {
    registerCommands(&baseCommands);
    registerCommands(&steamCommands);
}
void (*initCommands)(void) = init;
void initParameterSingleton() { new LocalParameters; }
