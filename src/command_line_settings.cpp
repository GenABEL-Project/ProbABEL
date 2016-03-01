/*
 * command_line_settings.cpp
 *
 *  Created on: Apr 1, 2012
 *      Author: mkooyman
*
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


#include <getopt.h>
#include <string>
#include <iostream>
#include "usage.h"
#include "command_line_settings.h"
#include "eigen_mematrix.h"

// config.h and fvlib/FileVector.h are included for the upper case variables
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "fvlib/FileVector.h"

using std::cout;
using std::cerr;
using std::endl;

string cmdvars::getStrGenfilename() const
{
    return str_genfilename;
}

int cmdvars::getAllcov() const
{
    return allcov;
}

string cmdvars::getChrom() const
{
    return chrom;
}

char* cmdvars::getGenfilename() const
{
    return genfilename;
}

int cmdvars::getInteraction() const
{
    return interaction;
}

char* cmdvars::getInverseFilename() const
{
    return inverse_filename;
}

bool cmdvars::isIscox() const
{
    return iscox;
}

int cmdvars::getIsFvf() const
{
    return isFVF;
}

char* cmdvars::getMapfilename() const
{
    return mapfilename;
}

char* cmdvars::getMlinfofilename() const
{
    return mlinfofilename;
}


int cmdvars::getNgpreds() const
{
    return ngpreds;
}

int cmdvars::getNohead() const
{
    return nohead;
}

int cmdvars::getNoutcomes() const
{
    return noutcomes;
}

int cmdvars::getNpeople() const
{
    return npeople;
}

string cmdvars::getOutfilename() const
{
    return outfilename;
}
//TODO(unknown) This function is not used. Remove in near future
//char* cmdvars::getProgramName() const
//{
//    return program_name;
//}

int cmdvars::getRobust() const
{
    return robust;
}

int cmdvars::getScore() const
{
    return score;
}

string cmdvars::getSep() const
{
    return sep;
}

int cmdvars::getSkipd() const
{
    return skipd;
}

char* cmdvars::getPhefilename() const
{
    return phefilename;
}


/**
 * Process the command line arguments and save them in a cmdvars object.
 *
 * @param argc Number of command line arguments
 * @param argv Values of the command line arguments
 */
void cmdvars::set_variables(int argc, char * argv[])
{
    int next_option;
    const char * const short_options = "p:i:d:m:n:c:o:s:t:g:a:relhb:v:u";
    // b - interaction parameter
    // ADD --fv FLAG (FILEVECTOR), IN WHICH CASE USE ALTERNATIVE
    // CONSTRUCTOR FOR GENDATA
    const struct option long_options[] =
    {
    { "pheno", 1, NULL, 'p' },
    { "info", 1, NULL, 'i' },
    { "dose", 1, NULL, 'd' },
    { "map", 1, NULL, 'm' },
    { "nids", 1, NULL, 'n' },
    { "chrom", 1, NULL, 'c' },
    { "out", 1, NULL, 'o' },
    { "skipd", 1, NULL, 's' },
    { "ntraits", 1, NULL, 't' },
    { "ngpreds", 1, NULL, 'g' },
    { "separat", 1, NULL, 'a' },
    { "score", 0, NULL, 'r' },
    { "no-head", 0, NULL, 'e' },
    { "allcov", 0, NULL, 'l' },
    { "help", 0, NULL, 'h' },
    { "interaction", 1, NULL, 'b' },
    { "mmscore", 1, NULL, 'v' },
    { "robust", 0, NULL, 'u' },
    { NULL, 0, NULL, 0 } };
    program_name = argv[0];

    do
    {
        next_option = getopt_long(argc, argv, short_options, long_options,
                NULL);

        switch (next_option)
        {
        case 'h':
            print_help(program_name, 0);
            break;
        case 'p':
            phefilename = optarg;
            neco[0] = 1;
            break;
        case 'i':
            mlinfofilename = optarg;
            neco[1] = 1;
            break;
        case 'd':
            genfilename = optarg;
            neco[2] = 1;
            break;
        case 'm':
            mapfilename = optarg;
            break;
        case 'n':
            npeople = atoi(optarg);
            break;
        case 'c':
            chrom = optarg;
            break;
        case 'o':
            outfilename = optarg;
            break;
        case 's':
            skipd = atoi(optarg);
            break;
        case 't':
            noutcomes = atoi(optarg);
            break;
        case 'g':
            ngpreds = atoi(optarg);
            break;
        case 'a':
            if (std::string(optarg) == std::string("\\t"))
            {
                sep = '\t';
            }
            else
            {
                sep = optarg;
            }
            break;
        case 'e':
            nohead = 1;
            break;
        case 'r':
            score = 1;
            break;
        case 'l':
            allcov = 1;
            break;
        case 'b':
            interaction = atoi(optarg);
            break;
        case 'v':
            inverse_filename = optarg;
            break;
        case 'u':
            robust = 1;
            break;

        case '?':
            print_usage(program_name, 2);
        case -1:
            break;
        default:
            abort();
        }  // end of switch
    } while (next_option != -1);  // end of while
}  // end of function


void cmdvars::printinfo()
{
    print_version();

    if (neco[0] != 1 || neco[1] != 1 || neco[2] != 1)
    {
        if (neco[0] != 1)
        {
            cerr << "Error: Missing required phenotype file (-p/--pheno option)"
                 << endl;
        }
        if (neco[1] != 1)
        {
            cerr << "Error: Missing required info file (-i/--info option)"
                 << endl;
        }
        if (neco[2] != 1)
        {
            cerr << "Error: Missing required genotype file (-d/--dose option)"
                 << endl;
        }
        cerr << endl;

        cout << "One or more required command line options "
             << "appear to be missing."
             << endl
             << "Run " << program_name
             << " --help for more information on the available options\n";
        exit(3);
    }


    if (score)
    {
        cout << "option --score suppressed from v 0.1-6\n";
        exit(1);
    }

    str_genfilename = genfilename;
    if (str_genfilename.find(FILEVECTOR_INDEX_FILE_SUFFIX) != string::npos
            || str_genfilename.find(FILEVECTOR_DATA_FILE_SUFFIX)
                    != string::npos)
        isFVF = 1;

    cout << "Options in effect:\n";
    cout << "\t --pheno       = "      << phefilename << endl;
    cout << "\t --info        = "      << mlinfofilename << endl;
    cout << "\t --dose        = "      << genfilename << endl;
    if (isFVF)
        cout << "\t             (using FVF data)" << endl;
    cout << "\t --ntraits     = "      << noutcomes << endl;
    cout << "\t --ngpreds     = "      << ngpreds << endl;
    cout << "\t --interaction = "      << interaction << endl;

    if (inverse_filename != NULL)
        cout << "\t --mmscore = "      << inverse_filename << endl;
    else
        cout << "\t --mmscore = not in output" << endl;

    if (mapfilename != NULL)
        cout << "\t --map     = "      <<  mapfilename << endl;
    else
        cout << "\t --map     = not in output" << endl;
    if (npeople > 0)
        cout << "\t --nids    = "      << npeople << endl;
    else
        cout << "\t --nids    = estimated from data" << endl;
    if (chrom != "-1")
        cout << "\t --chrom   = "      << chrom << endl;
    else
        cout << "\t --chrom   = not in output\n";
    if (outfilename.compare("") != 0)
        cout << "\t --out     = "      << outfilename << endl;
    else
        cout << "\t --out     = "      << "regression.out.txt" << endl;
    cout << "\t --skipd   = "          << skipd << endl;
    cout << "\t --separat = \""        << sep << "\"" << endl;
    if (score)
        cout << "\t --score   = ON"    << endl;
    else
        cout << "\t --score   = OFF"   << endl;
    if (nohead)
        cout << "\t --nohead  = ON"    << endl;
    else
        cout << "\t --nohead  = OFF"   << endl;
    if (allcov)
        cout << "\t --allcov  = ON"    << endl;
    else
        cout << "\t --allcov  = OFF"   << endl;
    if (robust)
        cout << "\t --robust  = ON"    << endl;
    else
        cout << "\t --robust  = OFF"   << endl;

    if (ngpreds != 1 && ngpreds != 2)
    {
        cerr << "\n\n--ngpreds should be 1 for MLDOSE or 2 for MLPROB"
             << endl;
        exit(1);
    }

    if (outfilename.compare("") == 0)
    {
        outfilename = string("regression");
    }
#if COXPH
    if (score)
    {
        cerr << "\n\nOption --score is implemented for "
             << "linear and logistic models only\n"
             << endl;
        exit(1);
    }

    if (inverse_filename != NULL)
    {
        cerr << "ERROR: mmscore is forbidden for cox regression"
             << endl;
        exit(1);
    }

    if (robust)
    {
        cerr << "ERROR: robust standard errors not implemented "
             << "for Cox regression"
             << endl;
        exit(1);
    }
#endif
    cout << endl;
}
