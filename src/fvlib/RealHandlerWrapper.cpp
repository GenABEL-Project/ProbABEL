#include "FileVector.h"
#include "RealHandlerWrapper.h"

// returns true on success
bool RealHandlerWrapper::open(const string &iFileName, bool iReadOnly) {
    fileName = iFileName;
    if (useCount > 0) {
        useCount++;
        return true;
    }

    if (iReadOnly) {
        stream.open(fileName.c_str(), ios::in | ios::binary);
    } else {
        stream.open(fileName.c_str(), ios::out | ios::in | ios::binary);
    }

    readOnly = iReadOnly;

    useCount = 1;

    return !!(stream);
}


void RealHandlerWrapper::close(){
    if (useCount > 1) {
        useCount--;
    } else if (useCount == 1) {
        useCount = 0;
        stream.close();
    }
}


void RealHandlerWrapper::blockWriteOrRead(unsigned long iLength, char* data,
                                          bool writeAction) {
    ::blockWriteOrRead(stream, iLength, data, writeAction);
}

void RealHandlerWrapper::fseek(unsigned long pos){
    stream.seekg(pos, ios::beg);
    stream.seekp(pos, ios::beg);
}

void RealHandlerWrapper::flush(){
    stream.flush();
}
