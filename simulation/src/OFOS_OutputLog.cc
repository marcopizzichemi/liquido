#include "OFOS_OutputLog.h" 

std::ofstream* OFOS_OutputLog::logfile = nullptr;
std::ofstream* OFOS_OutputLog::geom_logfile = nullptr;
std::stringstream OFOS_OutputLog::log_cache;
std::stringstream OFOS_OutputLog::ls_cache;
std::stringstream OFOS_OutputLog::geom_cache;
