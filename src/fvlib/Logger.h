#ifndef __LOGGER__
#define __LOGGER__

#if COMPILE_WITH_R
#include <R.h>
#endif

#include <ostream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string.h>

using namespace std;

#define MESSAGE_LEVEL 0
#define DEBUG_LEVEL 1
#define ERROR_LEVEL 2

class ErrorExit {};
class Logger {
 public:
//    enum ErrorLevel {Message, Debug, Error};
    int errorLevel;
    bool enabled;

Logger(int iErrorLevel) : errorLevel(iErrorLevel), enabled(true) {}
Logger(int iErrorLevel, bool iEnabled) : errorLevel(iErrorLevel),
        enabled(iEnabled) {}
    Logger &operator << (const char* s){
        sendString(s);
        return *this;
    }
    Logger &operator << (int x){ stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (unsigned int x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (long x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (void*x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (unsigned long x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (float x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (double x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (string x) { stringstream ss; string s; ss << x; ss >> s; sendString(s); return *this; }
    Logger &operator << (ostream&(*f)(ostream&)){ // endl support
        sendString("\n");
        return *this;
    }
    Logger &operator << (ErrorExit&) {
#ifdef R_R_H
        throw 1;
#else
        exit(EXIT_FAILURE);
#endif
        return *this;
    }

 private:
    void sendString(string s){
        if (!enabled)
            return;
#ifdef R_R_H
        Rprintf("%s",s.c_str());
#else
        if (errorLevel == ERROR_LEVEL) {
            cerr << s;
        } else {
            cout << s;
        }
#endif
    }
};

extern ErrorExit errorExit;
extern Logger inf;
extern Logger dbg;
extern Logger msg;
extern Logger errorLog;
extern Logger testDbg;
extern Logger deepDbg;
extern Logger fmDbg;
extern Logger wrapperLog;

#endif
