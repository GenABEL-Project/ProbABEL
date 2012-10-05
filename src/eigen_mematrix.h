#ifndef __EIGEN_MEMATRIX_H__
#define __EIGEN_MEMATRIX_H__
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

template<class DT> class mematrix
{
public:
    int nrow;
    int ncol;
    int nelements;
    Matrix<DT, Dynamic, Dynamic,RowMajor> data;

    mematrix()
    {
        nrow = ncol = nelements = 0;
//        std::cout << "[DEBUG]   creating object" << std::endl;

        data.resize(1, 1);
    }
    mematrix(int nr, int nc);
    mematrix(const mematrix &M);
    ~mematrix()
    {
//        if (nelements > 0)
//            delete data;
    }

    mematrix & operator=(const mematrix &M);
     DT &  operator[]( int i);
//    mematrix operator+(DT toadd);
    mematrix operator+(const mematrix &M);
    mematrix operator-(DT toadd);
    mematrix operator-(const mematrix &M);
    mematrix operator*(DT toadd);
    mematrix operator*(const mematrix &M);
    mematrix operator*(const mematrix *M);


    void reinit(int nr, int nc);

    unsigned int getnrow(void)
    {
        return nrow;
    }
    unsigned int getncol(void)
    {
        return ncol;
    }
    DT get(int nr, int nc);
    void put(DT value, int nr, int nc);
    DT column_mean(int nc);
    DT column_sum(int nc);
   void print(void);


};

//mematrix transpose(mematrix M);
//mematrix invert(mematrix M);

#endif
