#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <memory.h>

using namespace std;

#include "FilteredMatrix.h"
#include "CastUtils.h"
#include "frutil.h"

unsigned long FilteredMatrix::getCacheSizeInMb() {
        return nestedMatrix->getCacheSizeInMb();
}


void FilteredMatrix::setCacheSizeInMb( unsigned long cachesizeMb ) {
    nestedMatrix->setCacheSizeInMb(cachesizeMb);
}


void FilteredMatrix::setUpdateNamesOnWrite(bool bUpdate) {
    nestedMatrix->setUpdateNamesOnWrite(bUpdate);
}


void FilteredMatrix::writeVariableName(unsigned long varIdx, FixedChar name) {
    nestedMatrix->writeVariableName(filteredToRealRowIdx[varIdx], name);
}


void FilteredMatrix::writeObservationName(unsigned long obsIdx, FixedChar name) {
    nestedMatrix->writeObservationName(filteredToRealColIdx[obsIdx], name);
}


FixedChar FilteredMatrix::readVariableName(unsigned long varIdx) {
    return nestedMatrix->readVariableName(filteredToRealRowIdx[varIdx]);
}


FixedChar FilteredMatrix::readObservationName(unsigned long obsIdx) {
    return nestedMatrix->readObservationName(filteredToRealColIdx[obsIdx]);
}


void FilteredMatrix::readVariable(unsigned long varIdx, void * outvec) {
    unsigned long i;
    fmDbg << "readVariable(" << varIdx << "), numObservations="
          << getNumObservations() <<";" << endl;
    char* data =
        new (nothrow) char[getElementSize() * nestedMatrix->getNumObservations()];

    nestedMatrix->readVariable(this->filteredToRealRowIdx[varIdx], data);

    for (i = 0; i < this->filteredToRealColIdx.size(); i++){
        memcpy(&((char*)outvec)[i * getElementSize()],
               &data[this->filteredToRealColIdx[i] * getElementSize()],
               getElementSize());
    }
    delete [] data;
}


void FilteredMatrix::readObservation(unsigned long obsIdx, void * outvec) {
    unsigned long i;
    fmDbg << "readObservation(" << obsIdx << ");" << endl;
    for (i = 0; i < getNumVariables(); i++){
        readElement(i, obsIdx, (char*)outvec + i * getElementSize());
    }
}


void FilteredMatrix::writeObservation(unsigned long obsIdx, void * invec) {
    unsigned long i;
    for (i = 0; i < getNumObservations(); i++){
        writeElement(i, obsIdx, (char*)invec + i * getElementSize());
    }
}


void FilteredMatrix::writeVariable(unsigned long varIdx, void *datavec) {
    unsigned long i;
    fmDbg << "FilteredMatrix.writeVariable(" << varIdx << ")" << endl;
    double p = (double)getNumObservations() / nestedMatrix->getNumObservations();

    if (p > WRITE_SPEED_PROPORTION) {
        char *ptr = new char[getElementSize() *
                             nestedMatrix->getNumObservations()];
        // no filter
        if (getNumObservations() != nestedMatrix->getNumObservations()) {
            nestedMatrix->readVariable(this->filteredToRealRowIdx[varIdx], ptr);
        }
        for (i = 0; i < getNumObservations(); i++){
            memcpy(&ptr[getElementSize() * this->filteredToRealColIdx[i]],
                   &((char*)datavec)[getElementSize() * i],
            getElementSize());
        }

        nestedMatrix->writeVariable(this->filteredToRealRowIdx[varIdx], ptr);
        delete[] ptr;
    }else {
        for (i = 0; i < getNumObservations(); i++){
            this->writeElement(varIdx, i, (char*)datavec + i * getElementSize());
        }
    }
}


void FilteredMatrix::readElement(unsigned long varIdx, unsigned long obsIdx,
                                 void * out) {
    fmDbg << "FilteredMatrix::readElement(" << varIdx << "," << obsIdx<<") = ";
    nestedMatrix->readElement(filteredToRealRowIdx[varIdx],
                              filteredToRealColIdx[obsIdx], out);
    fmDbg << bufToString(getElementType(), (char*)out, string("NAN")) << endl;
}


void FilteredMatrix::writeElement(unsigned long varIdx, unsigned long obsIdx,
                                  void* data) {
    fmDbg << "FilteredMatrix.writeElement (" << varIdx << "," << obsIdx << ")"
          << endl;
    nestedMatrix->writeElement(filteredToRealRowIdx[varIdx],
                               filteredToRealColIdx[obsIdx], data);
}


string FilteredMatrix::getFileName(){
    return nestedMatrix->getFileName();
}


unsigned long FilteredMatrix::getNumVariables(){
   return filteredToRealRowIdx.size();
}


unsigned long FilteredMatrix::getNumObservations() {
   return filteredToRealColIdx.size();
}


void FilteredMatrix::saveAs(string newFilename) {
    nestedMatrix->saveAs(newFilename, this->filteredToRealRowIdx.size(),
                         this->filteredToRealColIdx.size(),
                         &this->filteredToRealRowIdx[0],
                         &this->filteredToRealColIdx[0]);
}


void FilteredMatrix::saveVariablesAs(string newFilename, unsigned long nvars,
                                     unsigned long *varIndexes) {
    vector<unsigned long> recodedColIndexes;
    vector<unsigned long> recodedRowIndexes;

    unsigned long *obsIndexes = new unsigned long[this->getNumObservations()];

    unsigned long i;
    for (i = 0; i < this->getNumObservations(); i++){
        obsIndexes[i]=i;
    }

    filterIdxList(obsIndexes, this->getNumObservations(), recodedColIndexes,
                  filteredToRealColIdx);
    filterIdxList(varIndexes, nvars, recodedRowIndexes, filteredToRealRowIdx);
    delete[] obsIndexes;
}


void FilteredMatrix::saveObservationsAs(string newFilename, unsigned long nobss,
                                        unsigned long * obsIndexes) {
    vector<unsigned long> recodedColIndexes;
    vector<unsigned long> recodedRowIndexes;

    unsigned long *varIndexes = new unsigned long[this->getNumVariables()];

    unsigned long i;
    for (i = 0; i < this->getNumObservations(); i++){
        varIndexes[i]=i;
    }

    filterIdxList(obsIndexes, nobss, recodedColIndexes, filteredToRealColIdx);
    filterIdxList(varIndexes, getNumVariables(), recodedRowIndexes,
                  filteredToRealRowIdx);
    delete obsIndexes;
}


void FilteredMatrix::saveAs(string newFilename, unsigned long nvars,
                            unsigned long nobss, unsigned long *varIndexes,
                            unsigned long *obsIndexes) {
    vector<unsigned long> recodedColIndexes;
    vector<unsigned long> recodedRowIndexes;
    filterIdxList(obsIndexes, nobss, recodedColIndexes, filteredToRealColIdx);
    filterIdxList(varIndexes, nvars, recodedRowIndexes, filteredToRealRowIdx);
    nestedMatrix->saveAs(newFilename, nvars, nobss, &recodedRowIndexes[0],
                         &recodedColIndexes[0]);
}


void FilteredMatrix::saveAsText(string newFilename, bool showVarNames,
                                bool showObsNames, string nanString) {
    nestedMatrix->saveAsText(newFilename, showVarNames, showObsNames, nanString);
}


short unsigned FilteredMatrix::getElementSize() {
    return nestedMatrix->getElementSize();
}


short unsigned FilteredMatrix::getElementType() {
    return nestedMatrix->getElementType();
}


void FilteredMatrix::addVariable(void * invec, string varname) {
    errorLog << "FilteredMatrix doesn't support addVariable."
             << endl << errorExit;
}


void FilteredMatrix::cacheAllNames(bool doCache) {
    nestedMatrix->cacheAllNames(doCache);
}


AbstractMatrix* FilteredMatrix::castToAbstractMatrix(){
    return this;
}


bool FilteredMatrix::setReadOnly(bool iReadOnly){
    return nestedMatrix->setReadOnly(iReadOnly);
}
