#ifndef __MEMATRIX_H__
#define __MEMATRIX_H__
#include <iostream>
using namespace std;

template<class DT> class mematrix {
public:
    int nrow;
    int ncol;
    int nelements;
    DT * data;

    mematrix() {
        nrow = ncol = nelements = 0;
        data = NULL;
    }

    mematrix(int nr, int nc);
    mematrix(const mematrix &M);
    ~mematrix() {
        if (nelements > 0)
            delete[] data;
    }

    mematrix & operator=(const mematrix &M);
    DT & operator[](int i);
    mematrix operator+(DT toadd);
    mematrix operator+(mematrix &M);
    mematrix operator-(DT toadd);
    mematrix operator-(mematrix &M);
    mematrix operator*(DT toadd);
    mematrix operator*(mematrix &M);

    void reinit(int nr, int nc);

    unsigned int getnrow(void) {
        return nrow;
    }
    unsigned int getncol(void) {
        return ncol;
    }
    DT get(int nr, int nc);
    void put(DT value, int nr, int nc);
    DT column_mean(int nc);
    void print(void);
    void delete_column(int delcol);
    void delete_row(int delrow);

};

//	mematrix transpose(mematrix M);
//	mematrix invert(mematrix M);

#endif
