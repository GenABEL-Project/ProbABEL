/*
 *
 * Copyright (C) 2009--2014 Various members of the GenABEL team. See
 * the SVN commit logs for more details.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301, USA.
 *
 */


#ifndef __MEMATRIX_H__
#define __MEMATRIX_H__
#include <iostream>
using namespace std;

template<class DT> class mematrix
{
 public:
    int nrow;
    int ncol;
    int nelements;
    DT * data;

    mematrix()
    {
        nrow = ncol = nelements = 0;
        data = NULL;
    }
    mematrix(int nr, int nc);
    mematrix(const mematrix &M);
    ~mematrix()
    {
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
    mematrix operator*(mematrix *M);

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
    void print(void);
    void delete_column(const int delcol);
    void delete_row(const int delrow);

};

//	mematrix transpose(mematrix M);
//	mematrix invert(mematrix M);

#endif
