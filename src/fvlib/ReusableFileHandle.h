#ifndef __ReusableFileHandle__
#define __ReusableFileHandle__

#include "RealHandlerWrapper.h"

#include <string>
#include <map>

using namespace std;

class ReusableFileHandle {
 private:
    bool isOk;
    unsigned long curPos;
    string fileName;
    bool readOnly;

    RealHandlerWrapper *realHandlerWrapper;

    static map<string, RealHandlerWrapper*> openHandles;

 public:
ReusableFileHandle(RealHandlerWrapper *iRealHandlerWrapper, bool iIsOk,
                   const string& iFileName, bool iReadOnly ):
    isOk(iIsOk), curPos(0),
        fileName(iFileName), readOnly(iReadOnly),
        realHandlerWrapper(iRealHandlerWrapper){
        }

ReusableFileHandle(const ReusableFileHandle &rfh) :
    isOk(rfh.isOk), curPos(rfh.curPos),fileName(rfh.fileName),
        readOnly(rfh.readOnly),    realHandlerWrapper(rfh.realHandlerWrapper)  {
    }

    void operator = (const ReusableFileHandle &rfh){
        isOk = rfh.isOk;
        curPos = rfh.curPos;
        fileName = rfh.fileName;
        readOnly = rfh.readOnly;
        realHandlerWrapper = rfh.realHandlerWrapper;
    }

ReusableFileHandle() : isOk(false),curPos(0),
        fileName(), readOnly(false)  {
    }

    operator bool(){
        return isOk;
    }

    void blockWriteOrRead(unsigned long length, char* data, bool writeAction);
    void fseek(unsigned long pos);
    void flush();

    static ReusableFileHandle getHandle(string fileName, bool readOnly);
    static int testGetNumHandles();
    void close();
};

#endif
