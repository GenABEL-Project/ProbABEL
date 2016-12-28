/**
 * \file   invsigma.h
 * \author mkooyman
 *
 * \brief Contains the definition of the InvSigma class.
 */
/*
 * Copyright (C) 2009--2016 Various members of the GenABEL team. See
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


#ifndef DATA_H_
#define DATA_H_
#include <string>

extern bool is_interaction_excluded;
//TODO(unknown) This function is not used. Remove in near future
//unsigned int Nmeasured(char * fname, int nphenocols, int npeople);
#include "phedata.h"
#include "gendata.h"


class InvSigma {
 private:
    unsigned int npeople;       /* number of people */
    std::string filename;
    mematrix<double> matrix;    /* file is stored here */

 public:
    InvSigma(const char * filename_, const phedata& phe);
    mematrix<double> & get_matrix();
    ~InvSigma();
};

#endif // DATA_H_
