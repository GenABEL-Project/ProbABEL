#ifndef __RealHandlerWrapper__
#define __RealHandlerWrapper__

class ReusableFileHandle;

#include <string>
#include <map>
#include <fstream>
#include <iostream>

using namespace std;

class RealHandlerWrapper {
 private:
    int useCount;
    fstream stream;
    string fileName;

    bool readOnly;

 public:
RealHandlerWrapper(): useCount(0) {}

    void blockWriteOrRead(unsigned long length, char* data, bool writeAction);
    void fseek(unsigned long pos);
    void flush();
    bool open(const string &fileName, bool readOnly);
    void close();

    int getUseCount() {return useCount;}
};

#endif
