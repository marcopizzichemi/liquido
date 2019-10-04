#ifndef OFOS_OutputLog_h
#define OFOS_OutputLog_h

#include <fstream>
#include <iostream>
#include <sstream>

class OFOS_OutputLog {
    public:
        static std::ofstream *logfile;
        static std::ofstream *geom_logfile;
        static std::stringstream log_cache;
        static std::stringstream ls_cache;
        static std::stringstream geom_cache;

};

#endif
