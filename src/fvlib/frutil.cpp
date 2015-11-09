#include <sys/stat.h>
#include <string>
#include <vector>

#include "frutil.h"
#include "const.h"
#include "FileVector.h"

FileHeader get_file_type(char * filename) {
    FileHeader out;
    ifstream myfile(filename, ios::binary | ios::in);
    if (!myfile) {
        errorLog << "can not open file for reading" << endl << errorExit;
    }
    myfile.read((char*)&out, sizeof(out));
    return(out);
}


string extract_base_file_name(string filename)
{
    unsigned int idxExtPos = filename.find(FILEVECTOR_INDEX_FILE_SUFFIX);
    unsigned int dataExtPos = filename.find(FILEVECTOR_DATA_FILE_SUFFIX);
    if (idxExtPos == filename.size() - FILEVECTOR_INDEX_FILE_SUFFIX.size() )
    {
        return filename.substr(0, idxExtPos);
    }
    else if (dataExtPos == filename.size() - FILEVECTOR_DATA_FILE_SUFFIX.size() )
    {
        return filename.substr(0, dataExtPos );
    }
    else
    {
        return filename;
    }
}


unsigned short calcDataSize(unsigned short int type){
    unsigned short desize;
    switch (type) {
    case UNSIGNED_SHORT_INT:
        desize = sizeof(unsigned short int);
        break;
    case SHORT_INT:
        desize = sizeof(short int);
        break;
    case UNSIGNED_INT:
        desize = sizeof(unsigned int);
        break;
    case INT:
        desize = sizeof(int);
        break;
    case FLOAT:
        desize = sizeof(float);
        break;
    case DOUBLE:
        desize = sizeof(double);
        break;
    case SIGNED_CHAR:
        desize = sizeof(char);
        break;
    case UNSIGNED_CHAR:
        desize = sizeof(unsigned char);
        break;
    default:
        desize = 0;
        errorLog << "file contains data of unknown type " << type
                 << endl << errorExit;
    }
    return desize;
}


void initializeEmptyFile(string filename, unsigned long numVariables,
                         unsigned long nobservations, unsigned short type,
                         bool override)
{
    dbg << "Initializing empty file '" << filename
        << "', type " << type << "." << endl;

    string indexFilename = filename + FILEVECTOR_INDEX_FILE_SUFFIX;
    string dataFilename = filename + FILEVECTOR_DATA_FILE_SUFFIX;

    FileHeader fileHeader;
    unsigned long desize = calcDataSize(type);
    fileHeader.type = type;
    fileHeader.numVariables = numVariables;
    fileHeader.numObservations = nobservations;
    fileHeader.nelements = numVariables * nobservations;
    fileHeader.bytesPerRecord = desize;
    fileHeader.bitsPerRecord = desize*8;

    bool bHeaderOrDataExists = headerOrDataExists( filename );

    if (override && bHeaderOrDataExists) {
        dbg << "Deleting existing file" << indexFilename << endl;
        unlink(indexFilename.c_str());
        unlink(dataFilename.c_str());
    }

    if (!override && bHeaderOrDataExists) {
        errorLog << "File '" << filename << "' already exists."
                 << endl << errorExit;
    }

    ofstream indexFile(indexFilename.c_str(), ios::binary | ios::out);
    ofstream dataFile(dataFilename.c_str(), ios::binary | ios::out);

    deepDbg << "Writing FileVector header." << endl;
    fileHeader.print();

    indexFile.seekp(0,ios::beg);
    indexFile.write((char*)&fileHeader, sizeof(fileHeader));

    deepDbg << "Writing " << nobservations << " observations." << endl;
    FixedChar name;
    for (unsigned long i = 0; i < nobservations; i++) {
        sprintf(name.name, "%lu", i+1);
        indexFile.seekp(sizeof(FileHeader) + i * sizeof(FixedChar), ios::beg);
        indexFile.write((char*)&name.name, sizeof(name.name));
    }

    deepDbg << "Writing " << numVariables << " variables." << endl;

    for (unsigned long j = 0; j < numVariables; j++) {
        sprintf(name.name, "%lu", j+1);
        indexFile.seekp(sizeof(FileHeader) +
                        (j + fileHeader.numObservations) * sizeof(FixedChar),
                        ios::beg);
        indexFile.write((char*)&name.name, sizeof(name.name));
    }
    indexFile.close();

    deepDbg << "Writing data file." << endl;

    unsigned long estimated_data_size =
        (unsigned long) fileHeader.bytesPerRecord *
        (unsigned long) fileHeader.numVariables *
        (unsigned long) fileHeader.numObservations;
    //unsigned long i;
    dataFile.seekp(estimated_data_size-1, ios::beg);
    char c = 0;
    dataFile.write(&c, 1);
    deepDbg << "Closing data file" << endl;
    dataFile.close();
    dbg << "File '" << filename << "' initialized." << endl;
}


bool headerOrDataExists(string fileName) {
    return file_exists(fileName + FILEVECTOR_INDEX_FILE_SUFFIX) ||
        file_exists(fileName + FILEVECTOR_INDEX_FILE_SUFFIX);
}


bool file_exists(string fileName) {
    struct stat buf;
    int i = stat( fileName.c_str(), &buf );
    /* File found */
    if (i == 0)
    {
        return true;
    }
    return false;
}


void tokenize(const string& str, vector<string>& tokens,
              const string& delimiters) {
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}


#define PART_SIZE INT_MAX
void blockWriteOrRead(fstream& file, unsigned long length,
                      char* data, bool writeAction){
    unsigned long i;
    unsigned long numParts = length/PART_SIZE;
    for(i = 0; i <= numParts; i++){
        unsigned long subLength =
            (numParts > 0 && i < numParts ) ? PART_SIZE:(length % PART_SIZE);
        if (writeAction){
            file.write(data+i*PART_SIZE, subLength);
        } else {
            file.read(data+i*PART_SIZE, subLength);
        }
    }
}
