#ifndef __FilteredMatrix__
#define __FilteredMatrix__

#include <map>

//todo add comments

using namespace std;

#include "AbstractMatrix.h"
/*
 * "Cols" are observations and "Rows" are variables
 */

class FilteredMatrix : public AbstractMatrix {
    AbstractMatrix *nestedMatrix;

    vector<unsigned long> filteredToRealColIdx;
    vector<unsigned long> filteredToRealRowIdx;

    void filterIdxList(unsigned long *iIndexes, unsigned long numIndexes,
                       vector<unsigned long> &oIndexes,
                       vector<unsigned long> &filter) {
        oIndexes.reserve(numIndexes);

        unsigned long i;
        for (i = 0; i < numIndexes; i++) {
            oIndexes.push_back(filter[iIndexes[i]]);
        }
    }

 public:
    // makes this matrix not filter any cells
    void setNoFiltering(){
        unsigned long i;
        filteredToRealRowIdx.reserve(nestedMatrix->getNumVariables());
        for (i = 0; i < nestedMatrix->getNumVariables(); i++) {
            filteredToRealRowIdx.push_back(i);
        }

        filteredToRealColIdx.reserve(nestedMatrix->getNumObservations());
        for (i = 0; i < nestedMatrix->getNumObservations(); i++) {
            filteredToRealColIdx.push_back(i);
        }
    }

    // set filter for Filterematrix
    void setFilteredArea(vector<unsigned long> &rowMask,
                         vector<unsigned long> &colMask){
        fmDbg << "setFilteredArea()" << endl;
        this->filteredToRealRowIdx = vector<unsigned long>(rowMask);
        this->filteredToRealColIdx = vector<unsigned long>(colMask);
    }

    /**
     * Constructs FilteredMatrix wrapper from AbstractMatrix object
     **/
FilteredMatrix(AbstractMatrix &matrix) : nestedMatrix(&matrix) {
        dbg << "Constructing FilteredMatrix from AbstractMatrix, ptr = "
            << (void *)this << endl;
        setNoFiltering();
        getWarningIsShown() = false;
    }

    string getFileName();

    /**
     * Returns number of variables
     **/
    unsigned long getNumVariables();

    /**
     * Returns number of observations
     **/
    unsigned long getNumObservations();

    /**
     * If argument is true, DatABEL reads all variables/observations'
     * names from files into memory.
     * Calling this function with false causes name cache freeing
     **/
    void cacheAllNames(bool);

    /**
     * Returns pointer to AbstractMatrix that is wrapped by current object
     **/
    AbstractMatrix* getNestedMatrix() {return nestedMatrix;}

    /**
     * Saves current DatABEL object to another file
     **/
    void saveAs(string newFilename);

    /**
     * Saves selected variables to another file
     **/
    void saveVariablesAs(string newFilename, unsigned long nvars,
                         unsigned long * varindexes);

    /**
     * Saves selected observations to another file
     **/
    void saveObservationsAs(string newFilename, unsigned long nobss,
                            unsigned long * obsindexes);

    /**
     * Saves selected observations to another file
     **/
    void saveAs(string newFilename, unsigned long nvars, unsigned long nobss,
                unsigned long * varindexes, unsigned long * obsindexes);

    /**
     * Saves current object as text
     **/
    void saveAsText(string newFilename, bool saveVarNames, bool saveObsNames,
                    string nanString);

    /**
     * Get single observation
     **/
    void readObservation(unsigned long obsIdx, void * outvec);

    /**
     * Write single observation
     **/
    void writeObservation(unsigned long obsIdx, void * invec);

    /**
     * Set variable name
     **/
    void writeVariableName(unsigned long varIdx, FixedChar newname);  // todo loooong future -- control that name is unique

    /**
     * Set observation name
     **/
    void writeObservationName(unsigned long obsIdx, FixedChar newname);  //todo loooong future -- control that name is unique!

    /**
     * Get cache size in Mbs
     **/
    unsigned long getCacheSizeInMb();

    /**
     * Set cache size in Mbs
     **/
    void setCacheSizeInMb( unsigned long cachesizeMb );

    /**
     * Get observation name
     **/
    FixedChar readObservationName(unsigned long obsIdx);

    /**
     * Get variable name
     **/
    FixedChar readVariableName(unsigned long varIdx);

    /**
     * If bUpdate is false and var/obs names are cached, changing var/obs name doesn't cause disk write.
     * Invoking this method with true forces update on each var/obs change.
     **/
    void setUpdateNamesOnWrite(bool bUpdate);

    /**
     * Returns number of bytes per element
     **/
    short unsigned getElementSize();

    /**
     * Returns elements type
     **/
    short unsigned getElementType();

    /**
     * Read a single variable
     **/
    void readVariable(unsigned long varIdx, void * outvec);

    /**
     * Read a single element
     **/
    void readElement(unsigned long varIdx, unsigned long obsIdx, void * elem);

    /**
     * Set a single variable
     **/
    void writeVariable(unsigned long varIdx, void * datavec);

    /**
     * Set a single element
     **/
    void writeElement(unsigned long varIdx, unsigned long obsIdx, void * data);

    /**
     * Returns ptr to current object as AbstractMatrix
     **/
    virtual AbstractMatrix* castToAbstractMatrix();

    /**
     * Enables readonly mode. Only one instance of Filtered matrix can
     * access data in write mode.
     **/
    virtual bool setReadOnly(bool readOnly);

 private:
    void addVariable(void * invec, string varname);
};

#endif
