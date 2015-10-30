#ifndef __FRUTIL__
#define __FRUTIL__

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "const.h"
#include "Logger.h"
#include "CastUtils.h"

using namespace std;

class FixedChar
{
 public:
    FixedChar(){
        memset(name,0xab,NAMELENGTH);
    };
    FixedChar(string s){
        if (s.length() >= NAMELENGTH){
            errorLog << "Overflow of FixedChar (length of name > NAMELENGTH ("
                     << NAMELENGTH << "): " << s.c_str() << "." << endl;

        }

        strncpy(name, s.c_str(), NAMELENGTH-1);
        name[NAMELENGTH-1] = '\0';
    };
    char name[NAMELENGTH];
};


class FileHeader
{
 public:
    unsigned short int type;
    // should change that to long!!!
    unsigned int nelements;
    unsigned int numObservations;
    unsigned int numVariables;
    unsigned int bytesPerRecord;
    unsigned int bitsPerRecord;
    unsigned int namelength;
    unsigned int reserved[RESERVEDSPACE];

    FileHeader () {
        type = 0;
        nelements = 0;
        numObservations = 0;
        bitsPerRecord = 0;
        bytesPerRecord = 0;
        numVariables = 0;
        namelength = NAMELENGTH;
        for (int i = 0; i < RESERVEDSPACE; i++) reserved[i] = 0;
    }
    ~FileHeader() {}

    void print() {
        dbg << "type = " << type << "("<< dataTypeToString(type) << ")" << endl;
        dbg << "nelements = " << nelements << endl;
        dbg << "numObservations = " << numObservations << endl;
        dbg << "numVariables = " << numVariables << ";" << endl;
        dbg << "bytesPerRecord = " << bytesPerRecord << ";" << endl;
        dbg << "bitsPerRecord = " << bitsPerRecord << ";" << endl;
    }
};


FileHeader get_file_type(char * filename);

void initializeEmptyFile(string filename, unsigned long numVariables,
                         unsigned long nobservations, unsigned short int type,
                         bool override);

string extract_base_file_name(string filename);
bool file_exists(string fileName);
bool headerOrDataExists(string fileName);
unsigned short calcDataSize(unsigned short int type);
void tokenize(const string& str, vector<string>& tokens,
              const string& delimiters = " \t");
void blockWriteOrRead(fstream& file, unsigned long length, char* data,
                      bool writeAction);
#endif
