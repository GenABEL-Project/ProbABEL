/*
 * gendata.cpp
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


#include <string>
#include <errno.h>
#include <limits>
#include "gendata.h"
#include "fvlib/FileVector.h"
#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif
#include "utilities.h"


void gendata::mldose_line_to_matrix(int k, const char *all_numbers,
                                    int amount_of_numbers){
    int j = 0;
    // Check if not a null pointer
    if (!*all_numbers){
        perror("Error while reading genetic data (expected pointer to char but found a null pointer)");
                       exit(EXIT_FAILURE);
    }

    while (j < amount_of_numbers)
    {
        double result = 0;
        // Skip whitespace
        while (*all_numbers == ' ')
        {
            all_numbers++;
        }

        // check NaN (right now checks only first character)
        // TODO: make catching of NaN more rigid
        if (*all_numbers == 'N')
        {
            result = std::numeric_limits<double>::quiet_NaN();
            // Skip other characters of NaN
            while ((*all_numbers == 'a') | (*all_numbers == 'N'))
            {
                all_numbers++;
            }
        }
        else
        {
            int sign = 0;
            // set sign to -1 if negative: multiply by sign just before return
            if (*all_numbers == '-')
            {
                all_numbers++;
                sign = -1;
            }
            // Read digits before dot
            while (*all_numbers <= '9' && *all_numbers >= '0')
            {
                result = result * 10 + (*all_numbers++ - '0');
            }
            // Read digit after dot
            if (*all_numbers == '.')
            {
                double decimal_counter = 1.0;
                all_numbers++;
                while (*all_numbers <= '9' && *all_numbers >= '0')
                {
                    decimal_counter *= 0.1;
                    result += (*all_numbers++ - '0') * decimal_counter;
                }
            }
            // Correct for negative number
            if (sign == -1)
            {
                result = sign * result;
            }
        }
        G.put(result, k, j);
        j++;
    }
}


void gendata::get_var(int var, double * data)
{
    // Read the genetic data for SNP 'var' and store in the array 'data'

    if (DAG == NULL)            // Read from text file
    {
        for (int i = 0; i < G.nrow; i++)
        {
            data[i] = G.get(i, var);
        }
    }
    else if (DAG != NULL)       // Read from fv file
    {
        // cout << "Data Type: " << dataTypeToString(DAG->getElementType())
        //      << endl;
        double *tmpdata = new double[DAG->getNumObservations()];
        DAG->readVariableAs((unsigned long int) var, tmpdata);

        unsigned int j = 0;
        for (unsigned int i = 0; i < DAG->getNumObservations(); i++)
        {
            if (!DAGmask[i])
            {
                // A dirty trick to get rid of conversion
                // errors. Instead of casting float data to double we
                // convert the data to string and then do strtod()
                char tmpstr[1048576];
                snprintf (tmpstr, sizeof(tmpstr), "%f", tmpdata[i]);

                double val;
                char *endptr;
                errno = 0;      // To distinguish success/failure
                                // after strtod()
                val = strtod(tmpstr, &endptr);

                if ((errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALL))
                    || (errno != 0 && val == 0)) {
                    perror("Error while reading genetic data (strtod)");
                    exit(EXIT_FAILURE);
                }

                if (endptr == tmpstr) {
                    cerr << "No digits were found while reading genetic data"
                         << " (individual " << i + 1
                         << ", position " << var + 1 << ")"
                         << endl;
                    exit(EXIT_FAILURE);
                }
                /* If we got here, strtod() successfully parsed a number */
                data[j++] = val;
            }
        }
        delete[] tmpdata;
    }
    else
    {
        report_error("cannot get gendata");
    }
}


gendata::gendata() : nsnps(0), nids(0), ngpreds(0), DAG(NULL), DAGmask(NULL)
{
}


void gendata::re_gendata(string filename, unsigned int insnps,
                         unsigned int ingpreds, unsigned int npeople,
                         unsigned int nmeasured,
                         unsigned short int * allmeasured,
                         std::string * idnames)
{
    nsnps = insnps;
    ngpreds = ingpreds;
    DAG = new FileVector(filename, 128, true);
    DAGmask = new unsigned short int[DAG->getNumObservations()];
    if (DAG->getNumObservations() != npeople)
        report_error("dimension of fvf-data and phenotype data do not match\n");

    if (DAG->getNumVariables() != insnps * ingpreds)
        report_error("dimension of fvf-data and mlinfo data do not match\n");

    long int j = -1;

    for (unsigned int i = 0; i < npeople; i++)
    {
        if (allmeasured[i] == 0)
        {
            DAGmask[i] = 1;
        }
        else
        {
            DAGmask[i] = 0;
            j++;
        }
        string DAGobsname = DAG->readObservationName(i).name;

        if (DAGobsname.find("->") != string::npos)
        {
            DAGobsname = DAGobsname.substr(DAGobsname.find("->") + 2);
        }

        // if (allmeasured[i] && idnames[j] != DAGobsname)
        //  std::cerr << "names do not match for observation at phenofile "
        //            << "line (phe/geno) " << i+1 << "/+1 ("
        //            << idnames[i].c_str() << "/"
        //            << DAGobsname.c_str() << ")\n";
        // fix thanks to Vadym Pinchuk
        if (allmeasured[i] && idnames[j] != DAGobsname)
        {
            report_error(
                "names do not match for observation at phenofile line (phe/geno) %i/+1 (%s/%s)\n",
                i + 1, idnames[j].c_str(), DAGobsname.c_str());
        }
    }
    nids = j + 1;
    // std::cout << "in INI: " << nids << " " << npeople << "\n";
    if (nids != nmeasured)
        report_error("nids != mneasured (%i != %i)\n", nids, nmeasured);
}


void gendata::re_gendata(char * fname, unsigned int insnps,
                         unsigned int ingpreds, unsigned int npeople,
                         unsigned int nmeasured,
                         unsigned short int * allmeasured, int skipd,
                         std::string * idnames)
{
    nids    = nmeasured;
    nsnps   = insnps;
    ngpreds = ingpreds;
    DAG     = NULL;
    //	int nids_all = npeople;


    G.reinit(nids, (nsnps * ngpreds));

    std::ifstream infile;
    infile.open(fname);

    if (!infile)
    {
        std::cerr << "gendata: cannot open file " << fname << endl;
    }

    std::string tmpid, tmpstr;

    int k = 0;
    for (unsigned int i = 0; i < npeople; i++)
    {
        if (allmeasured[i] == 1)
        {
            if (skipd > 0)
            {
                // Read the genotype data and look for the signature
                // arrow of MaCH/minimac. If found only use the part
                // after the arrow as ID.
                infile >> tmpstr;
                size_t strpos = tmpstr.find("->");
                if (strpos != string::npos)
                {
                    tmpid = tmpstr.substr(strpos + 2, string::npos);
                }
                else
                {
                    tmpid = tmpstr;
                }
                if (tmpid != idnames[k])
                {
                    cerr << "phenotype file and dose or probability file "
                            << "did not match at line " << i + 2 << " ("
                            << tmpid << " != " << idnames[k] << ")" << endl;
                    infile.close();
                    exit(1);
                }
            }

            for (int j = 1; j < skipd; j++)
            {
                infile >> tmpstr;
            }

            std::string all_numbers;
            all_numbers.reserve(nsnps * ngpreds * 7);
            std::getline(infile, all_numbers);
            mldose_line_to_matrix(k, all_numbers.c_str(), nsnps * ngpreds);

            k++;
        }
        else
        {
            for (int j = 0; j < skipd; j++)
            {
                infile >> tmpstr;
            }
            for (unsigned int j = 0; j < (nsnps * ngpreds); j++)
            {
                infile >> tmpstr;
            }
        }
    }

    infile.close();

}

// HERE NEED A NEW CONSTRUCTOR BASED ON DATABELBASECPP OBJECT
gendata::~gendata()
{
    if (DAG != NULL)
    {
        delete DAG;
        delete[] DAGmask;
    }

    //		delete G;
}
