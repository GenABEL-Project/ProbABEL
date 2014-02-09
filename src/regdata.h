/*
 * regdata.h
 *
 *  Created on: Mar 29, 2012
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


#ifndef REGDATA_H_
#define REGDATA_H_

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif
#include "gendata.h"
#include "phedata.h"

class regdata {
 public:
    int nids;
    int ncov;
    int ngpreds;
    int noutcomes;
    bool is_interaction_excluded;
    unsigned short int * masked_data;
    unsigned int gcount;
    double freq;
    mematrix<double> X;
    mematrix<double> Y;
    regdata();
    regdata(const regdata &obj);
    regdata(phedata &phed, gendata &gend, const int snpnum,
            const bool ext_is_interaction_excluded);
    mematrix<double> extract_genotypes();
    void update_snp(gendata *gend, const int snpnum);
    void remove_snp_from_X();
    regdata get_unmasked_data();
    ~regdata();

 private:
};

#endif /* REGDATA_H_ */
