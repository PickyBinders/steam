#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "createteasubdb.sh.h"

int createteasubdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program = par.db3 + ".sh";
    FileUtil::writeFile(program, createteasubdb_sh, createteasubdb_sh_len);

    CommandCaller cmd;
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("LIST", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("CREATESUBDB1_PAR", par.createParameterString(par.createsubdb).c_str());
    par.dbIdMode = 0;
    cmd.addVariable("CREATESUBDB2_PAR", par.createParameterString(par.createsubdb).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    assert(false);
    return EXIT_FAILURE;
}
