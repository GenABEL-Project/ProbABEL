#include "Logger.h"


ErrorExit errorExit;

#define all true

Logger inf(MESSAGE_LEVEL,all);
Logger dbg(DEBUG_LEVEL,false);
Logger msg(MESSAGE_LEVEL,all);
Logger testDbg(DEBUG_LEVEL,all);
Logger deepDbg(DEBUG_LEVEL,false);
Logger fmDbg(DEBUG_LEVEL,false);
Logger wrapperLog(DEBUG_LEVEL,false);
Logger errorLog(ERROR_LEVEL,all);
