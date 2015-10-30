#ifndef __FRVECTOR__
#define __FRVECTOR__

#include <new>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>

#include "Logger.h"
#include "frutil.h"

#include "ReusableFileHandle.h"
#include "AbstractMatrix.h"

using namespace std;

#define DB_CREATE 1
#define DB_EXCL 2
#define DB_RDONLY 4


class FileVector: public AbstractMatrix {
 private:
    const string filename;
    string dataFilename;
    string indexFilename;
    ReusableFileHandle dataFile;
    ReusableFileHandle indexFile;
    FileHeader fileHeader;
    // row and column names
    FixedChar *variableNames;
    FixedChar *observationNames;
    // size of header (descriptives + var/obs names)
    unsigned long headerSize;
    // cache size (Mb) requested by user
    unsigned long cache_size_Mb;
    // cache size internal; these ones are exact and used internaly
    unsigned long cache_size_nvars;
    unsigned long cache_size_bytes;
    unsigned long max_buffer_size_bytes;
    // which variables are now in cache
    unsigned long cacheBegin;
    unsigned long cacheEnd;
    char * cacheBuffer;
    bool readOnly;
    bool updateNamesOnWrite;

 public:
    FileVector();
    ~FileVector();

FileVector(string iFilename, unsigned long cachesizeMb) : filename (iFilename) {
        readOnly = false;
        updateNamesOnWrite = false;
        cacheBuffer = 0;
        initialize(cachesizeMb);
        getWarningIsShown() = false;
    }

FileVector(string iFilename, unsigned long cachesizeMb, bool iReadOnly) :
    filename(iFilename), readOnly(iReadOnly) {
        updateNamesOnWrite = false;
        cacheBuffer = 0;
        initialize(cachesizeMb);
        getWarningIsShown() = false;
    }

FileVector(char *iFilename, unsigned long cachesizeMb) :
    filename(string(iFilename)){
        updateNamesOnWrite = false;
        readOnly = false;
        cacheBuffer = 0;
        initialize(cachesizeMb);
        getWarningIsShown() = false;
    }

FileVector(char *iFilename, unsigned long cachesizeMb, bool iReadOnly) :
    filename(string(iFilename)), readOnly(iReadOnly) {
        updateNamesOnWrite = false;
        cacheBuffer = 0;
        initialize(cachesizeMb);
        getWarningIsShown() = false;
    }
    // for testing purposes
    void getPrivateCacheData(unsigned long* cacheSizeNVars,
                             unsigned long *inCachFrom, unsigned long *cacheEnd);

    // these ones are the actual used to initialize and free up
    void initialize(unsigned long cachesizeMb);
    void deInitialize();
    // this one updates cache
    void updateCache(unsigned long from_var);

    // gives element number from varIdx & obsIdx
    unsigned long nrnc_to_nelem(unsigned long varIdx, unsigned long obsIdx);

    // getting and setting var/col names
    void writeVariableName(unsigned long varIdx, FixedChar name);
    void writeObservationName(unsigned long obsIdx, FixedChar name);

    virtual string getFileName();
    virtual unsigned long getNumVariables();
    virtual unsigned long getNumObservations();

    FixedChar readVariableName(unsigned long varIdx);
    FixedChar readObservationName(unsigned long obsIdx);

    // USER FUNCTIONS
    // can read single variable
    void readVariable(unsigned long varIdx, void * outvec);

    void addVariable(void * invec, string varname);

    // should only be used for reading single random elements!
    void readElement(unsigned long varIdx, unsigned long obsIdx, void * data);
    void writeVariable(unsigned long varIdx, void * datavec);
    void writeElement(unsigned long varIdx, unsigned long obsIdx, void * data);

    void readObservation(unsigned long obsIdx, void * outvec);
    void writeObservation(unsigned long obsIdx, void * outvec);
    void cacheAllNames(bool);
    void readNames();

    void saveAs( string newFilename );
    void saveVariablesAs(string newFilename, unsigned long nvars,
                         unsigned long * varindexes);
    void saveObservationsAs(string newFilename, unsigned long nobss,
                            unsigned long * obsindexes);
    void saveAs(string newFilename, unsigned long nvars, unsigned long nobss,
                unsigned long * varindexes, unsigned long * obsindexes);
    void saveAsText(string newFilename, bool saveVarNames, bool saveObsNames,
                    string nanString);

    unsigned long getCacheSizeInMb();
    void setCacheSizeInMb( unsigned long cachesizeMb );

    virtual short unsigned getElementSize();
    virtual short unsigned getElementType();
    void saveIndexFile();

    virtual void setUpdateNamesOnWrite(bool bUpdate);
    virtual AbstractMatrix* castToAbstractMatrix();
    virtual bool setReadOnly(bool readOnly);

    // FOR FUTURE:
    // very slow one!
    //	DT * readObservation(unsigned long obsIdx);
    // should only be used for reading single random elements!
    //	DT readElement(unsigned long nelment);
 private :
    void copyVariable(char * to, char * from, int n, unsigned long * indexes );
    void calcCachePos(unsigned long newPos, unsigned long &cacheBegin,
                      unsigned long &cacheEnd);
};

//global variables
const string FILEVECTOR_DATA_FILE_SUFFIX = ".fvd";
const string FILEVECTOR_INDEX_FILE_SUFFIX = ".fvi";

#endif
