#include <sys/stat.h>
#include <string.h>
#include <iostream>

using namespace std;

#include "FileVector.h"
#include "frutil.h"


void FileVector::saveIndexFile() {
    if (readOnly)
        return;

    indexFile.fseek(0);
    indexFile.blockWriteOrRead(sizeof(fileHeader), (char*)&fileHeader, true);
    indexFile.fseek(sizeof(fileHeader));

    if (observationNames && variableNames) {
        indexFile.blockWriteOrRead(sizeof(FixedChar) * fileHeader.numObservations,
                                   (char*)observationNames, true);
        indexFile.fseek(sizeof(fileHeader) + sizeof(FixedChar) *
                        fileHeader.numObservations);
        indexFile.blockWriteOrRead(sizeof(FixedChar) * fileHeader.numVariables,
                                   (char*)variableNames, true);
    }
}


void FileVector::deInitialize(){
    saveIndexFile();
    delete [] cacheBuffer;
    cacheBuffer = 0;
    delete [] observationNames;
    observationNames = 0;
    delete [] variableNames;
    variableNames = 0;
    indexFile.close();
    dataFile.close();
    AbstractMatrix::closeForWriting(filename);
}


FileVector::~FileVector() {
    deInitialize();
}


void FileVector::initialize(unsigned long cachesizeMb) {
    //Rprintf("initialize, open file %s\n",filename.c_str());
    dbg << "Opening FileVector '" << filename.c_str() <<"'."<< endl;

    if (!readOnly) {
        AbstractMatrix::checkOpenForWriting(filename);
    }

    indexFilename = extract_base_file_name(filename) +
        FILEVECTOR_INDEX_FILE_SUFFIX;
    dataFilename = extract_base_file_name(filename) +
        FILEVECTOR_DATA_FILE_SUFFIX;

    if (!file_exists(indexFilename)) {
        errorLog << "Index file not exists: " <<  indexFilename
                 << endl << errorExit;
    }

    dataFilename = extract_base_file_name(filename) + FILEVECTOR_DATA_FILE_SUFFIX;
    if (!file_exists(dataFilename))
        errorLog << "Data file not exists: " <<  dataFilename.c_str()
                 << endl <<errorExit;

    struct stat data_filestatus;
    stat(dataFilename.c_str(), &data_filestatus);

    struct stat index_filestatus;
    stat(indexFilename.c_str(), &index_filestatus);

    indexFile = ReusableFileHandle::getHandle(indexFilename, readOnly);

    if (!indexFile) {
        errorLog << "Opening file " << indexFilename
                 << " for write & read failed\n" << errorExit;
    }

    dataFile = ReusableFileHandle::getHandle(dataFilename, readOnly);

    if (!dataFile) {
        errorLog << "Opening file " << dataFilename
                 << " for write & read failed\n" << errorExit;
    }

    indexFile.blockWriteOrRead(sizeof(fileHeader), (char*)&fileHeader, false);
    if (!indexFile){
        errorLog << "Failed to read datainfo from file:"
                 << indexFilename << endl;
    }

    // some integrity checks
    if (getElementSize() != fileHeader.bytesPerRecord) {
        errorLog << "System data type size (" << getElementSize();
        errorLog <<") and file data type size ("
                 << fileHeader.bytesPerRecord <<") do not match.\n";
    }

    //!!! nelements should actually be long to ensure !!!
    if (fileHeader.nelements !=
        (fileHeader.numObservations * fileHeader.numVariables)) {
        errorLog << "Number of variables (" << fileHeader.numVariables;
        errorLog << ") and observations (" << fileHeader.numObservations
                 << ") do not multiply to nelements";
        errorLog << "(" << fileHeader.nelements << ") (file integrity issue?)\n";
        errorLog << errorExit;
    }

    if ((fileHeader.bytesPerRecord != (fileHeader.bitsPerRecord/8)) ||
        ((fileHeader.bitsPerRecord % 8) != 0) || (fileHeader.bitsPerRecord < 8)) {
        errorLog << "Size in bytes/bits do not match or bit-size of char !=8 or ";
        errorLog << "non-byte recods (file integrity issue?)" << errorExit;
    }

    unsigned long indexSize = sizeof(fileHeader) + sizeof(FixedChar) *
        (fileHeader.numVariables + fileHeader.numObservations);
    if(indexSize != (unsigned long) index_filestatus.st_size) {
        errorLog << "Index file " << indexFilename << " size("
                 << (int) index_filestatus.st_size
                 << ") differs from the expected(";
        errorLog << indexSize <<")" << endl << errorExit;
    }

    // temp fix because nelements is not yet long ... !!!
    //	unsigned long estimated_size = data_type.bytes_per_record*data_type.nelements + headerSize;
    unsigned long estimated_size =
        (unsigned long) fileHeader.bytesPerRecord *
        (unsigned long) fileHeader.numVariables *
        (unsigned long) fileHeader.numObservations;

    if (estimated_size != (unsigned long) data_filestatus.st_size) {
        errorLog << "Data file size (" << (int) data_filestatus.st_size;
        errorLog << ") differs from the expected (" << estimated_size << ")"
                 << endl << " [";
        errorLog << fileHeader.numVariables << ","
                 << fileHeader.numObservations<<"]" << endl;
        errorLog << errorExit;
    }

    variableNames = 0;
    observationNames = 0;

    setCacheSizeInMb(cachesizeMb);
    updateCache(0);
    dbg << "Filevector " << filename << " opened." << endl;
}


// read in variable and observation names
void FileVector::readNames() {
    if (variableNames) delete[] variableNames;
    if (observationNames) delete[] observationNames;

    variableNames = new (nothrow) FixedChar [fileHeader.numVariables];
    if (!variableNames)
        errorLog << "can not get RAM for variable names" << errorExit;

    observationNames = new (nothrow) FixedChar [fileHeader.numObservations];
    if (!observationNames)
        errorLog << "can not get RAM for observation names" << errorExit;

    indexFile.fseek(sizeof(fileHeader));
    for (unsigned long i = 0; i < fileHeader.numObservations; i++)
    {
        indexFile.blockWriteOrRead(sizeof(FixedChar),
                                   (char*)(observationNames+i), false);
    }

    for (unsigned long i = 0; i < fileHeader.numVariables; i++)
        indexFile.blockWriteOrRead(sizeof(FixedChar),
                                   (char*)(variableNames+i), false);
}

unsigned long FileVector::getCacheSizeInMb() {
    return cache_size_Mb;
}

void FileVector::setCacheSizeInMb( unsigned long cachesizeMb ) {
    // figure out cache size
    cache_size_Mb = cachesizeMb;
    cache_size_nvars = (unsigned long) 1024 * 1024 * cache_size_Mb /
        (fileHeader.numObservations * fileHeader.bytesPerRecord);

    if (cache_size_nvars < 1) {
        cache_size_Mb = (unsigned long) ceil(
            (float) fileHeader.numObservations * fileHeader.bytesPerRecord /
            (1024.*1024.)
            );
        cache_size_nvars = 1;
    } else if (cache_size_nvars > fileHeader.numVariables) {
        cache_size_Mb = (unsigned long) ceil(
            (float) fileHeader.numVariables * fileHeader.numObservations *
            fileHeader.bytesPerRecord / (1024. * 1024.)
            );
        cache_size_nvars = fileHeader.numVariables;
    }
    cache_size_bytes = cache_size_nvars * fileHeader.bytesPerRecord *
        fileHeader.numObservations * sizeof(char);

    //free previously allocated memory
    if(cacheBuffer !=0 )
        delete[] cacheBuffer;
    // get memory for the cache
    cacheBuffer = new (nothrow) char [cache_size_bytes];
    if (!cacheBuffer)
        errorLog << "failed to get memory for cache" << endl << errorExit;

    max_buffer_size_bytes = INT_MAX;

    //don't read cache after resizing,
    //it will be updated on next read operation from desired position
    cacheBegin = 1;
    cacheEnd = 0;
}


void FileVector::calcCachePos(unsigned long newCenterPos,
                              unsigned long &cacheBeginRef,
                              unsigned long &cacheEndRef){
    if (cache_size_nvars == getNumVariables()) {
        cacheBeginRef = 0;
        cacheEndRef = getNumVariables();
        return;
    }

    cacheBeginRef = newCenterPos - cache_size_nvars/2;
    cacheEndRef = cacheBeginRef + cache_size_nvars;

    //cacheBeginRef < 0 ?
    if (newCenterPos < cache_size_nvars / 2) {
        cacheBeginRef = 0;
        cacheEndRef = cacheBeginRef + cache_size_nvars;
        return;
    }

    if (cacheEndRef > getNumVariables()){
        cacheEndRef = getNumVariables();
        cacheBeginRef = cacheEndRef - cache_size_nvars;
    }
}


void FileVector::updateCache(unsigned long varIdx) {
    // first load ?
    if (cacheEnd == 0 && cacheBegin == 1){
        calcCachePos(varIdx, cacheBegin, cacheEnd);
        dataFile.fseek(cacheBegin);
        dbg << "First time cache load." << endl;
        dataFile.blockWriteOrRead(cache_size_bytes, cacheBuffer, false);
        if (!dataFile){
            errorLog << "Inner error reading file." << endl << errorExit;
        }
        return;
    }

    if (getNumObservations() == 0) {
        return;
    }

    unsigned long currentCenter = (cacheEnd + cacheBegin) / 2;
    unsigned long newCenter = varIdx;

    unsigned long diffBetweenCenters;

    if (currentCenter>newCenter){
        diffBetweenCenters = currentCenter - newCenter;
    } else {
        diffBetweenCenters = newCenter - currentCenter;
    }

    if (diffBetweenCenters < cache_size_nvars/4 ){
        // no need to update cache
        return;
    }

    unsigned long newCacheBegin;
    unsigned long newCacheEnd;

    calcCachePos(newCenter, newCacheBegin, newCacheEnd);

    // update cache ?
    if (newCacheBegin == cacheBegin) {
        return;
    }

    unsigned long oldPos;
    unsigned long newPos;
    unsigned long loadPos;
    unsigned long inMemLoadPos;
    unsigned long hddReadSize;
    unsigned long memMoveSize;

    if (newCacheBegin > cacheBegin){
        oldPos = newCacheBegin - cacheBegin;
        newPos = 0;
        loadPos = max(newCacheBegin, cacheEnd);
        inMemLoadPos = max(newCacheBegin, cacheEnd) - newCacheBegin;
        hddReadSize = min(newCacheBegin, cacheEnd ) - cacheBegin;
    } else {
        oldPos = 0;
        newPos = cacheBegin - newCacheBegin;
        loadPos = newCacheBegin;
        inMemLoadPos = 0;
        hddReadSize = min(newCacheEnd, cacheBegin) - newCacheBegin;
    }

    memMoveSize = cache_size_nvars - hddReadSize;

    if (memMoveSize>0){
        memmove(cacheBuffer + newPos * getElementSize() * getNumObservations(),
                cacheBuffer + oldPos * getElementSize() * getNumObservations(),
                memMoveSize * getElementSize() * getNumObservations());
    }

    dataFile.fseek(loadPos * getElementSize() * getNumObservations());
    dataFile.blockWriteOrRead(hddReadSize * getElementSize() *
                              getNumObservations(),
                              cacheBuffer + inMemLoadPos * getElementSize() *
                              getNumObservations(), false);
    if (!dataFile){
        errorLog << "Inner error reading file." << endl << errorExit;
    }

    cacheBegin = newCacheBegin;
    cacheEnd = newCacheEnd;
}


void FileVector::setUpdateNamesOnWrite(bool bUpdate) {
    updateNamesOnWrite = bUpdate;
}


void FileVector::cacheAllNames(bool doCache) {
    if (doCache) {
        if (variableNames || observationNames) {
            dbg << "FileVector.cacheAllNames(true) called while variable "
                << "names are already cached." << endl;
            return;
        }
        readNames();
    } else {
        if (variableNames) {
            delete[] variableNames;
            variableNames = 0;
        }
        if (observationNames) {
            delete[] observationNames;
            observationNames = 0;
        }
    }
}


void FileVector::writeVariableName(unsigned long varIdx, FixedChar name) {
    if (varIdx >= fileHeader.numVariables) {
        errorLog << "Trying to set name of obs out of range (" << varIdx
                 << ")\n\n" << endl << errorExit;
    }

    if ((updateNamesOnWrite || variableNames == 0) && !readOnly){
        indexFile.fseek(sizeof(fileHeader) + sizeof(FixedChar) *
                        (varIdx + fileHeader.numObservations));
        indexFile.blockWriteOrRead(sizeof(FixedChar), (char*)&name, true);
        indexFile.flush();
    }

    if (variableNames) {
        variableNames[varIdx] = name;
    }
}


void FileVector::writeObservationName(unsigned long obsIdx, FixedChar name) {
    if (obsIdx >= fileHeader.numObservations) {
        errorLog << "Trying to set name of vars out of range (" << obsIdx
                 << ")\n\n" << endl << errorExit;
    }

    if ((updateNamesOnWrite || observationNames == 0) && !readOnly){
        indexFile.fseek(sizeof(fileHeader) + sizeof(FixedChar) * (obsIdx));
        indexFile.blockWriteOrRead(sizeof(FixedChar), (char*)&name, true);
        indexFile.flush();
    }

    if (observationNames) {
        observationNames[obsIdx] = name;
    }
}


FixedChar FileVector::readVariableName(unsigned long varIdx) {
    if (varIdx >= fileHeader.numVariables) {
        errorLog << "trying to get name of var out of range" << errorExit;
    }

    if (!variableNames) {
        FixedChar ret;
        indexFile.fseek(sizeof(fileHeader) +
                        sizeof(FixedChar) * (varIdx+fileHeader.numObservations));
        indexFile.blockWriteOrRead(sizeof(FixedChar), (char*)&ret, false);
        return ret;
    }

    return variableNames[varIdx];
}


FixedChar FileVector::readObservationName(unsigned long obsIdx) {
    if (obsIdx >= fileHeader.numObservations) {
        errorLog << "trying to get name of obs out of range" << errorExit;
    }

    if (!observationNames) {
        FixedChar ret;
        indexFile.fseek(sizeof(fileHeader) + sizeof(FixedChar)*(obsIdx));
        indexFile.blockWriteOrRead(sizeof(FixedChar), (char*)&ret, false);
        return ret;
    }

    return observationNames[obsIdx];
}


// can read single variable
void FileVector::readVariable(unsigned long varIdx, void * outvec) {
    if (varIdx >= fileHeader.numVariables) {
        errorLog << "Variable number out of range (" << varIdx
                 << " >= " << fileHeader.numVariables << ")"
                 << endl << errorExit;
    }
    updateCache(varIdx);
    unsigned long offset = (varIdx - cacheBegin) *
        fileHeader.numObservations * getElementSize();
    memcpy(outvec, cacheBuffer+offset,
           getElementSize() * (fileHeader.numObservations));
}


void FileVector::readObservation(unsigned long obsIdx, void *outvec) {
    char * tmpdata = new (nothrow) char [getNumObservations()*getElementSize()];
    if (!tmpdata)
        errorLog << "readObservation: cannot allocate tmpdata" << errorExit;

    for (unsigned long int i = 0; i< getNumVariables(); i++)
    {
        readVariable(i, tmpdata);
        memcpy((char*)outvec + i * getElementSize(),
               tmpdata + getElementSize() * obsIdx,
               getElementSize());
    }
    delete[] tmpdata;
}


void FileVector::writeObservation(unsigned long obsIdx, void * invec) {
    if (readOnly) {
        errorLog << "Trying to write to the readonly file." << errorExit;
    }
    for (unsigned long int i = 0; i < getNumVariables(); i++)	{
        writeElement(i, obsIdx, (char*)invec + i * getElementSize() );
    }
}


// can write single variable
void FileVector::writeVariable(unsigned long varIdx, void * datavec) {
    if (readOnly) {
        errorLog << "Trying to write to the readonly file." << errorExit;
    }
    unsigned long pos = nrnc_to_nelem(varIdx, 0);
    dataFile.fseek(pos * getElementSize());
    dataFile.blockWriteOrRead(getElementSize() * fileHeader.numObservations,
                              (char*)datavec, true);
    dataFile.flush();
    if (!dataFile) {
        errorLog << "failed to write to data file\n" << errorExit;
    }

    if (varIdx >= cacheBegin && varIdx < cacheEnd)
    {
        unsigned long offset = (varIdx - cacheBegin) *
            fileHeader.numObservations * getElementSize();
        memcpy(cacheBuffer + offset,
               datavec,
               getElementSize() * fileHeader.numObservations);
    }
}


unsigned long FileVector::nrnc_to_nelem(unsigned long varIdx,
                                        unsigned long obsIdx) {
    if (varIdx >= fileHeader.numVariables) {
        errorLog << "Variable number out of bounds (" << varIdx
                 << " >= " <<  fileHeader.numVariables << ")"
                 << endl << errorExit;
    }

    if (obsIdx >= fileHeader.numObservations) {
        errorLog << "Observation number out of bounds (" << obsIdx
                 << " >= " <<  fileHeader.numVariables << ")"
                 << endl << errorExit;
    }
    return( varIdx * fileHeader.numObservations + obsIdx );
}


// should only be used for reading single random elements!
void FileVector::readElement(unsigned long varIdx, unsigned long obsIdx,
                             void* out) {
    unsigned long pos = nrnc_to_nelem(varIdx, obsIdx);
    deepDbg << "FileVector.readElement(" << varIdx << "," << obsIdx
            << "), pos = " << pos << ", ";
    dataFile.fseek(pos * getElementSize());
    dataFile.blockWriteOrRead(getElementSize(), (char*)out, false);
}


void FileVector::writeElement(unsigned long varIdx, unsigned long obsIdx,
                              void* data) {
    if (readOnly) {
        errorLog << "Trying to write to the readonly file." << errorExit;
    }
    deepDbg << "FileVector.writeElement(" << varIdx << "," << obsIdx << ");"
            << endl;
    unsigned long pos = nrnc_to_nelem(varIdx, obsIdx);
    dataFile.fseek(pos * getElementSize());
    dataFile.blockWriteOrRead(getElementSize(), (char*)data, true);
    dataFile.flush();

    if (varIdx >= cacheBegin && varIdx < cacheEnd) {
        unsigned long offset = (varIdx - cacheBegin) *
            fileHeader.numObservations * getElementSize() +
            obsIdx *getElementSize();
        memcpy(cacheBuffer+offset, data, getElementSize() );
    }
}


unsigned long FileVector::getNumVariables() {
    return fileHeader.numVariables;
}


unsigned long FileVector::getNumObservations() {
    return fileHeader.numObservations;
}


void FileVector::saveAs( string newFilename ) {
    initializeEmptyFile( (char *)newFilename.c_str(),
                         getNumVariables(), getNumObservations(),
                         fileHeader.type, true);
    FileVector *outdata = new FileVector( newFilename, 64 );  //todo which size for cache to use?

    // copy observation names from the first object
    for (unsigned long i = 0; i < getNumObservations(); i++)
        outdata->writeObservationName( i, readObservationName(i) );

    char * tmpvariable =
        new (nothrow) char[getNumObservations() * getElementSize()];

    if (!tmpvariable)
        errorLog << "can not allocate memory for tmpvariable"
                 << endl << endl << errorExit;

    for (unsigned long i = 0; i < getNumVariables(); i++)
    {
        //write var names
        outdata->writeVariableName(i, readVariableName(i));
        //write variables
        readVariable(i, tmpvariable);
        outdata->writeVariable(i, tmpvariable);
    }
    delete outdata;
    delete [] tmpvariable;
}


void FileVector::saveVariablesAs(string newFilename, unsigned long nvars,
                                 unsigned long *varindexes) {
    initializeEmptyFile( (char *)newFilename.c_str(), nvars,
                         getNumObservations(), fileHeader.type, true);
    FileVector outdata( newFilename, 64 );  //todo which size for cache to use?

    // copy observation names from the first object
    for (unsigned long i = 0; i < getNumObservations(); i++)
        outdata.writeObservationName(i, readObservationName(i) );

    char * tmpvariable =
        new (nothrow) char[getNumObservations() * getElementSize()];

    if (!tmpvariable)
        errorLog << "can not allocate memory for tmpvariable"
                 << endl << endl << errorExit;

    for (unsigned long i = 0; i < nvars; i++ )
    {
        unsigned long selected_index =  varindexes[i];
        //write var names
        outdata.writeVariableName(i, readVariableName(selected_index));
        //write variables
        readVariable(selected_index, tmpvariable);
        outdata.writeVariable(i, tmpvariable);
    }

    delete [] tmpvariable;
}


string FileVector::getFileName() {
    return filename;
}


void FileVector::saveObservationsAs( string newFilename, unsigned long nobss,
                                     unsigned long *obsindexes) {
    if (headerOrDataExists(newFilename))
    {
        errorLog << "File " << newFilename <<" already exists"
                 << endl << errorExit;
    }

    initializeEmptyFile( (char *)newFilename.c_str(), getNumVariables(),
                         nobss, fileHeader.type, true);
    FileVector outdata( newFilename, 64 );  //todo which size for cache to use?

    // copy observation names from the first object
    for (unsigned long i = 0; i < nobss; i++ )
        outdata.writeObservationName(i, readObservationName( obsindexes[i] ) );

    char * in_variable =
        new (nothrow) char[getNumObservations() * getElementSize()];

    if (!in_variable)
        errorLog << "can not allocate memory for tmpvariable"
                 << endl << endl << errorExit;

    char * out_variable = new (nothrow) char[nobss*getElementSize()];
    if (!out_variable)
        errorLog << "can not allocate memory for tmpvariable"
                 << endl << endl << errorExit;

    for (unsigned long i = 0; i < getNumVariables(); i++)
    {
        //write var names
        outdata.writeVariableName(i, readVariableName(i));
        //write variables
        readVariable(i, in_variable);
        copyVariable(out_variable, in_variable, nobss, obsindexes);
        outdata.writeVariable(i, out_variable);
    }

    delete [] in_variable;
    delete [] out_variable;
}


/*
 * copy elements from "from" array to "to" array, according to "n" and "indexes" parameters
 */
void FileVector::copyVariable(char* to, char* from, int n,
                              unsigned long * indexes ) {
    for (int j = 0; j < n; j++ ) {
        //copy only selected observations to out_variable  from in_variable
        unsigned long int read_offset = indexes[j] * getElementSize();
        if(read_offset + getElementSize() >
           getNumObservations() * getElementSize()) {
            errorLog << "When saving selected observations: index in obsindexes("
                     << indexes[j];
            errorLog << ") is out of range, source obsIdx is "
                     << getNumObservations() << endl;
            errorLog << errorExit;
        }
        memcpy(to + j*getElementSize(),from + read_offset,getElementSize());
    }
}


void FileVector::saveAs(string newFilename, unsigned long nvars,
                        unsigned long nobss, unsigned long *varindexes,
                        unsigned long *obsindexes) {
    if (headerOrDataExists(newFilename)) {
        errorLog << "File "<< newFilename <<" already exists."
                 << endl << errorExit;
    }
    initializeEmptyFile( (char *)newFilename.c_str(), nvars, nobss,
                         fileHeader.type, true);
    FileVector outdata(newFilename, 64);  //todo which size for cache to use?

    // copy observation names from the first object
    for (unsigned long i = 0; i < nobss; i++) {
        //cout << nobss << " " << i << "-" << obsindexes[i] << ";";
        outdata.writeObservationName(i, readObservationName(obsindexes[i]) );
    }

    char * out_variable = new (nothrow) char[nobss*getElementSize()];
    if (!out_variable)
        errorLog << "can not allocate memory for out_variable"
                 << endl << errorExit;

    char * in_variable =
        new (nothrow) char[getNumObservations()*getElementSize()];

    if (!in_variable)
        errorLog << "can not allocate memory for in_variable"
                 << endl << errorExit;

    for (unsigned long i= 0; i < nvars; i++) {
        unsigned long selected_index =  varindexes[i];
        //write var names
        outdata.writeVariableName(i, readVariableName(selected_index));
        //write variables
        readVariable(selected_index,in_variable);
        copyVariable(out_variable, in_variable, nobss, obsindexes);
        outdata.writeVariable(i,out_variable );
    }

    delete[] in_variable;
    delete[] out_variable;
}


void FileVector::saveAsText(string newFilename, bool saveVarNames,
                            bool saveObsNames, string nanString) {

    ofstream textfile(newFilename.c_str(), ios::out);

    // copy observation names from the first object
    if (saveObsNames) {
        for (unsigned long i = 0; i < getNumObservations(); i++ ) {
            FixedChar fc = readObservationName(i) ;
            textfile << fc.name << " ";
        }

        textfile << endl;
    }

    char * in_variable =
        new (nothrow) char[getNumObservations() * getElementSize()];

    if (!in_variable)
        errorLog << "can not allocate memory for in_variable"
                 << endl << endl << errorExit;

    for(unsigned long i = 0; i < getNumVariables(); i++) {
        dbg << "Writing var " << i << " of " << getNumVariables() << endl;
        //write var names
        FixedChar fc = readVariableName(i);
        if (saveVarNames) {
            textfile << fc.name << " ";
        }
        //write variables
        readVariable(i, in_variable);

        for (unsigned long j = 0; j < getNumObservations(); j++) {
            string s  = bufToString(getElementType(),
                                    &in_variable[j * getElementSize()],
                                    nanString);
            textfile << s << " ";
        }
        textfile << endl;
    }

    delete[] in_variable;
}


short unsigned FileVector::getElementSize() {
    return calcDataSize(fileHeader.type);
}


short unsigned FileVector::getElementType() {
    return fileHeader.type;
}


void FileVector::addVariable(void *invec, string varName) {
    deepDbg << "addVariable(" << varName << ")" << endl;
    if (readOnly) {
        errorLog << "Trying to write to the readonly file." << errorExit;
    }

    fileHeader.numVariables++;
    //recalculate
    fileHeader.nelements = fileHeader.numVariables * fileHeader.numObservations;

    FixedChar _fc_varname(varName);

    // are names loaded from disk ?
    if (variableNames && observationNames) {

        FixedChar * newVariablesNames =
            new (nothrow)FixedChar[fileHeader.numVariables];

        if (!newVariablesNames) {
            errorLog << "Can not allocate memory in addVariable()" << errorExit;
        }

        //reallocate greater array for var names
        memcpy(newVariablesNames, variableNames,
               sizeof(FixedChar)*(fileHeader.numVariables-1));
        newVariablesNames[fileHeader.numVariables - 1] = _fc_varname;
        FixedChar *oldvar_names = variableNames;
        variableNames = newVariablesNames;
        delete[] oldvar_names;

        if (updateNamesOnWrite) {
            indexFile.fseek(sizeof(fileHeader) +
                            sizeof(FixedChar) * (fileHeader.numVariables - 1 +
                                                 fileHeader.numObservations));
            indexFile.blockWriteOrRead(sizeof(FixedChar),
                                       (char*)&_fc_varname.name, true);
        }
    } else { // not loaded
        indexFile.fseek(sizeof(fileHeader) +
                        sizeof(FixedChar) * (fileHeader.numVariables - 1 +
                                             fileHeader.numObservations));
        indexFile.blockWriteOrRead(sizeof(FixedChar),
                                   (char*)&_fc_varname.name, true);
    }
    writeVariable(fileHeader.numVariables - 1, invec);
}


void FileVector::getPrivateCacheData(unsigned long* cacheSizeNVars,
                                     unsigned long *pCacheBegin,
                                     unsigned long *pCacheEnd ) {
    *cacheSizeNVars = cache_size_nvars;
    *pCacheBegin = cacheBegin;
    *pCacheEnd = cacheEnd;
}


bool FileVector::setReadOnly(bool iReadOnly){
    if (iReadOnly) {
        if (!this->readOnly) {
            deInitialize();
            this->readOnly = iReadOnly;
            initialize(this->cache_size_Mb);
        }
    } else {
        if (this->readOnly) {
            bool canOpen;
            {
                ofstream indexFileTest(indexFilename.c_str(),
                                       ios::out|ios::in|ios::binary);
                ofstream dataFileTest(dataFilename.c_str(),
                                      ios::out|ios::in|ios::binary);
                canOpen = indexFileTest.good() && dataFileTest.good();
            }

            if (canOpen) {
                deInitialize();
                this->readOnly = iReadOnly;
                initialize(this->cache_size_Mb);
            } else {
                errorLog << "Can't open " << filename << "for writing. " << endl;
                return false;
            }
        }
    }
    return true;
}


AbstractMatrix* FileVector::castToAbstractMatrix(){
    return this;
}
