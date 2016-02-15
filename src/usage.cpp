/*
 *
 * Copyright (C) 2009--2015 Various members of the GenABEL team. See
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


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "command_line_settings.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "eigen_mematrix.h"


using std::cout;
using std::endl;


void print_usage(char * program_name, int exit_code)
{
    cout << "\nUsage: " << program_name << " options" << endl;
    cout << "Options:" << endl;
    cout << "\t --pheno   : phenotype file name" << endl;
    cout << "\t --info    : information (e.g. MLINFO) file name" << endl;
    cout << "\t --dose    : predictor (e.g. MLDOSE/MLPROB) file name"
         << endl;
    cout << "\t --map     : [optional] map file name" << endl;
    cout << "\t --nids    : [optional] number of people to analyse" << endl;
    cout << "\t --chrom   : [optional] chromosome (to be passed to output)"
         << endl;
    cout << "\t --out     : [optional] output file name "
         << "(default is regression.out.txt)"
         << endl;
    cout << "\t --skipd   : [optional] how many columns to skip in the "
         << "predictor\n"
         << "\t             (dose/prob) file (default 2)"
         << endl;
#if COXPH
    cout << "\t --ntraits : [optional] how many traits are analysed (default 2)"
         << endl;
#else
    cout << "\t --ntraits : [optional] how many traits are analysed (default 1)"
         << endl;
#endif
    cout << "\t --ngpreds : [optional] how many predictor columns per marker\n"
         <<"\t              (default 1 = MLDOSE; else use 2 for MLPROB)"
         << endl;
    cout << "\t --separat : [optional] character to separate fields "
         << "(default is space)"
         << endl;
    cout << "\t --mafflip : swap/flip reference and effect allele based on MAF"
         << endl;
    cout << "\t --score   : use score test" << endl;
    cout << "\t --no-head : do not report header line" << endl;
    cout << "\t --allcov  : report estimates for all covariates "
         << "(large outputs!)"
         << endl;
    cout << "\t --interaction: Which covariate to use for interaction with "
         << "SNP analysis (default is no interaction, 0)"
         << endl;
    cout << "\t --interaction_only: like previous but without covariate acting"
         << " in interaction with SNP (default is no interaction, 0)"
         << endl;
#if !COXPH
    cout << "\t --mmscore : score test in samples of related individuals. "
         << "File with inverse of variance-covariance matrix (for palinear)"
         << " or inverse correlation (for palogist) as input parameter"
         << endl;
    cout << "\t --robust  : report robust (aka sandwich, aka Hubert-White) "
         << "standard errors"
         << endl;
#endif
    cout << "\t --help    : print help" << endl;
    exit(exit_code);
}


void print_version(void) {
    cout << PACKAGE
         << " v. " << PACKAGE_VERSION
         << "\n(C) Yurii Aulchenko, Lennart C. Karssen, Maarten Kooyman, "
         << "Maksim Struchalin, The GenABEL team, EMC Rotterdam\n\n";
    cout << "Using EIGEN version " << EIGEN_WORLD_VERSION
         << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION
         << " for matrix operations\n";
}


void print_help(char * program_name, int exit_code)
{
    print_version();
    print_usage(program_name, exit_code);
}
