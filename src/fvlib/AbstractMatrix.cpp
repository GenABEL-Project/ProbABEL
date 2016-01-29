#include "AbstractMatrix.h"

set<string> AbstractMatrix::fileNamesOpenForWriting;

void AbstractMatrix::checkOpenForWriting(const string fileName){
    deepDbg << "checkOpenForWriting(" << fileName << ")" << endl;
    if (AbstractMatrix::fileNamesOpenForWriting.find(fileName) !=
        fileNamesOpenForWriting.end()) {
        errorLog << "File " << fileName << " is already opened." <<  endl;
        throw 1;
    } else {
        AbstractMatrix::fileNamesOpenForWriting.insert(fileName);
    }
}


void AbstractMatrix::closeForWriting(const string fileName){
    fmDbg << "closeForWriting(" << fileName << ")" << endl;
    AbstractMatrix::fileNamesOpenForWriting.erase(fileName);
}
