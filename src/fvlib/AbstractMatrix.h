#ifndef __AbstractMatrix__
#define __AbstractMatrix__

#include <string>
#include <set>

using namespace std;

#include "frutil.h"
#include "CastUtils.h"

#define WRITE_SPEED_PROPORTION .01

// See filteredMatrix.h for detailed comments

class AbstractMatrix {
 public:
    virtual ~AbstractMatrix(){};

    template <class DT>
        void writeVariableAs(unsigned long varIdx, DT * outvec)
    {
        char* tmp = new (nothrow) char [getNumObservations()*getElementSize()];
        if (!tmp)
            errorLog << "writeVariableAs allocation error" << errorExit;
        for (unsigned long int i = 0; i < getNumObservations(); i++){
            performCast(&tmp[i * getElementSize()], outvec[i], getElementType(),
                        warningIsShown);
        }
        writeVariable(varIdx, tmp);
        delete[] tmp;
    }

    template <class DT>
        void addVariableAs(DT * outvec, string varname)
    {
        char* tmp = new (nothrow) char [getNumObservations()*getElementSize()];
        if (!tmp)
            errorLog << "add_variable_as allocation error" << errorExit;
        for (unsigned long int i = 0; i < getNumObservations(); i++){
            performCast(&tmp[i * getElementSize()], outvec[i], getElementType(),
                        warningIsShown);
        }
        addVariable (tmp,  varname);
        delete[] tmp;
    }

    template<class DT>
        void readVariableAs(unsigned long varIdx, DT * outvec)
    {
        char * tmp = new (nothrow) char[getNumObservations()*getElementSize()];
        readVariable(varIdx, tmp);
        for (unsigned long int i = 0; i < getNumObservations(); i++) {
            performCast(outvec[i], &tmp[i*getElementSize()], getElementType(),
                        warningIsShown);
        }
        delete[] tmp;
    }

    template<class DT>
        void readElementAs(unsigned long varNumber, unsigned long obsNumber,
                           DT & element){
        char *ret= new char [getElementSize()];
        readElement(varNumber, obsNumber, ret);
        performCast(element, ret, getElementType(), warningIsShown);
        delete [] ret;
    }

    template <class DT>
        void writeElementAs(unsigned long varNumber, unsigned long obsNumber,
                            DT& element){
        deepDbg << "AbstractMatrix.writeElementAs(" << varNumber << ","
                << obsNumber << "," << element <<")";
        deepDbg << "Alloc getElementSize() = " << getElementSize() << endl;
        char *ret = new char [getElementSize()];
        deepDbg << "Perform cast" << endl;
        performCast(ret, element, getElementType(), warningIsShown);
        writeElement(varNumber, obsNumber, ret);
        delete [] ret;
    }

    virtual string getFileName() = 0;

    virtual unsigned long getNumVariables() = 0;
    virtual unsigned long getNumObservations() = 0;

    virtual void saveAs(string newFilename) = 0;
    virtual void saveVariablesAs(string newFilename, unsigned long nvars,
                                 unsigned long * varindexes) = 0;
    virtual void saveObservationsAs(string newFilename, unsigned long nobss,
                                    unsigned long * obsindexes) = 0;

    virtual void saveAs(string newFilename, unsigned long nvars,
                        unsigned long nobss, unsigned long * varindexes,
                        unsigned long * obsindexes) = 0;
    virtual void saveAsText(string newFilename, bool saveVarNames,
                            bool saveObsNames, string nanString) = 0;

    virtual void readObservation(unsigned long obsIdx, void * outvec) = 0;
    virtual void writeObservation(unsigned long obsIdx, void * invec) = 0;

    virtual void writeVariableName(unsigned long varIdx, FixedChar newname) = 0;    // todo loooong future -- control that name is unique
    virtual void writeObservationName(unsigned long obsIdx, FixedChar newname)= 0;   //todo loooong future -- control that name is unique!

    virtual unsigned long getCacheSizeInMb() = 0;
    virtual void setCacheSizeInMb( unsigned long cachesizeMb ) = 0;

    virtual FixedChar readObservationName(unsigned long obsIdx) = 0;
    virtual FixedChar readVariableName(unsigned long varIdx) = 0;
    virtual void cacheAllNames(bool) = 0;

    virtual void setUpdateNamesOnWrite(bool bUpdate) = 0;
    virtual short unsigned getElementSize() = 0;
    virtual short unsigned getElementType() = 0;
    virtual void readVariable(unsigned long varIdx, void * outvec) = 0;
    virtual void readElement(unsigned long varIdx, unsigned long obsIdx,
                             void * elem) = 0;
    virtual void writeVariable(unsigned long varIdx, void * datavec) = 0;
    virtual void writeElement(unsigned long varIdx, unsigned long obsIdx,
                              void * data) = 0;
    virtual AbstractMatrix* castToAbstractMatrix() = 0;
    virtual bool setReadOnly(bool readOnly) = 0;

    static set<string> fileNamesOpenForWriting;
    static void checkOpenForWriting(const string fileName);
    static void closeForWriting(const string fileName);

    bool &getWarningIsShown(){ return warningIsShown;}
 private:

    // HIGH -- here I see the possibility to make these functions faster then "random" access functions
    // adds variable at the end = writeVariable with varIdx=NVARS?
    // todo loooong future -- control that name is unique!
    virtual void addVariable(void * invec, string varname) = 0;
    //    virtual void add_observation(void * invec, string obsname) = 0;
    // write single element
    // CURRENTLY CACHE IS NOT UPDATED!
    bool warningIsShown;
};

#endif
