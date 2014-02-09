/**
 * \file   data.h
 * \author mkooyman
 *
 * \brief Contains several classes we didn't put somewhere else yet
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


#ifndef DATA_H_
#define DATA_H_
#include <string>

extern bool is_interaction_excluded;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople);
#include "phedata.h"
#include "gendata.h"

/**
 * \brief Data from the mlinfo file.
 *
 */
class mlinfo {
 public:
    int nsnps;                  /**< Number of SNPs */
    std::string * name;         /**< Array of SNP names */
    std::string * A1;           /**< Array with the first allele */
    std::string * A2;           /**< Array with the second allele */
    double * Freq1;
    double * MAF;               /**< The minor allele frequency */
    double * Quality;           /**< The imputation quality metric */
    double * Rsq;               /**< The imputation \f$R^2\f$ */
    std::string * map;          /**< Array with the SNP positions */
    mlinfo()
    {
        Freq1 = NULL;
        MAF = NULL;
        Quality = NULL;
        Rsq = NULL;
        nsnps = 0;
        A1 = NULL;
        A2 = NULL;
        name = NULL;
        map = NULL;
    }
    mlinfo(char * filename, char * mapname);
    ~mlinfo();
};

class InvSigma {
 private:
    static const unsigned MAXIMUM_PEOPLE_AMOUNT = 1000000;
    unsigned int npeople;       /* number of people */
    std::string filename;
    mematrix<double> matrix;    /* file is stored here */

 public:
    InvSigma(const char * filename_, phedata * phe);
    mematrix<double> & get_matrix();
    ~InvSigma();
};

#endif /* DATA_H_ */
