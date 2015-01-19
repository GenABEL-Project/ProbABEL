/**
 * \file   mlinfo.h
 * \author The GenABEL team
 *
 * \brief Contains the class information for the mlinfo class.
 *
 * The mlinfo class contains the information read from the .info of
 * .mlinfo files.
 *
 */
/* Copyright (C) 2009--2015 Various members of the GenABEL team. See
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

#ifndef MLINFO_H_
#define MLINFO_H_
#include <string>
#include <cstdlib>

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


    // Constructors and destructors
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

#endif // MLINFO_H_
