/*
 * gendata.h
 *
 *  Created on: Mar 8, 2012
 *      Author: mkooyman
 *
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


#ifndef GENDATA_H_
#define GENDATA_H_
#include <string>
#include "fvlib/FileVector.h"

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#endif

class gendata {
 public:
    unsigned int nsnps;
    unsigned int nids;
    unsigned int ngpreds;
    gendata();

    void re_gendata(char * fname, unsigned int insnps, unsigned int ingpreds,
            unsigned int npeople, unsigned int nmeasured,
            unsigned short int * allmeasured, int skipd, std::string * idnames);

    void re_gendata(string filename, unsigned int insnps, unsigned int ingpreds,
            unsigned int npeople, unsigned int nmeasured,
            unsigned short int * allmeasured, std::string * idnames);

    void get_var(int var, double * data);

    ~gendata();

    // MAKE THAT PRIVATE, ACCESS THROUGH GET_SNP
    // ANOTHER PRIVATE OBJECT IS A POINTER TO DATABELBASECPP
    // UPDATE SNP, ALL REGRESSION METHODS: ACCOUNT FOR MISSING
 private:
    mematrix<double> G;
    AbstractMatrix * DAG;
    unsigned short int * DAGmask;
};

#endif /* GENDATA_H_ */
