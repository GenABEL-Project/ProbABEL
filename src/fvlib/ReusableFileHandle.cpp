#include <stdlib.h>

#include "RealHandlerWrapper.h"
#include "ReusableFileHandle.h"


map<string, RealHandlerWrapper*> ReusableFileHandle::openHandles;

int ReusableFileHandle::testGetNumHandles(){
    return ReusableFileHandle::openHandles.size();
}


ReusableFileHandle ReusableFileHandle::getHandle(string fileName, bool readOnly){
    string key = string((readOnly ? "R" : "*")) + fileName;

    if (ReusableFileHandle::openHandles.find(key) !=
        ReusableFileHandle::openHandles.end()) {
        RealHandlerWrapper *rhw = ReusableFileHandle::openHandles[key];

        rhw->open(fileName, readOnly);

        ReusableFileHandle ret(rhw, true, fileName, readOnly);

        return ret;
    } else {
        RealHandlerWrapper *newHandleWrapper = new RealHandlerWrapper();
        bool success = newHandleWrapper->open(fileName, readOnly);

        if (success) {
            ReusableFileHandle::openHandles[key] = newHandleWrapper;
        } else {
            delete newHandleWrapper;
            newHandleWrapper = 0;
        }

        ReusableFileHandle ret(newHandleWrapper, success, fileName, readOnly);
        return ret;
    }
}


void ReusableFileHandle::close() {
    string key = string((readOnly ? "R" : "*")) + fileName;

    if (ReusableFileHandle::openHandles.find(key) ==
        ReusableFileHandle::openHandles.end())
        return;

    RealHandlerWrapper *rhw = ReusableFileHandle::openHandles[key];

    rhw->close();

    if (rhw->getUseCount() == 0) {
        delete rhw;
        ReusableFileHandle::openHandles.erase(key);
    }
}


void ReusableFileHandle::blockWriteOrRead(unsigned long length, char* data,
                                          bool writeAction){
    realHandlerWrapper->fseek(curPos);
    realHandlerWrapper->blockWriteOrRead(length, data, writeAction);

    curPos += length;
}


void ReusableFileHandle::flush(){
    realHandlerWrapper->flush();
}


void ReusableFileHandle::fseek(unsigned long pos){
    curPos = pos;
}
