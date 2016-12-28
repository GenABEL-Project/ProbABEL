/**
 * \file phedata.cpp
 * \author Y.S. Aulchenko
 * \author L.C. Karssen
 * \author M. Kooyman
 * \author Maksim V. Struchalin
 *
 * \brief Contains the functions of the phdata class containing the
 * phenotype data.
 *
 *
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


#include <phedata.h>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cstdlib>

using std::cout;
using std::cerr;
using std::endl;


phedata::phedata(const char * fname, const int noutc, const int npeople,
                 const int interaction, const bool iscox)
{
    setphedata(fname, noutc, npeople, interaction, iscox);
}


/**
 * Read phenotype data from file.
 *
 * @param fname Name of the file containing phenotype data
 * @param noutc Number of outcomes/phenotypes in the phenotype file
 * (see #noutcomes).
 * @param npeople Number of people to use in the analysis. If set to
 * 0, then all individuals in the phenotype file will be used. If > 0
 * (i.e. set using the \--nids command line option, see
 * cmdvars::set_variables()) only the first npeople will be used.
 * @param interaction Column specifying which phenotype/covariate is
 * selected to interact with the SNP (default: 0, i.e. no interaction).
 * @param iscox Are we running a Cox PH regression?
 */
void phedata::setphedata(const char * fname, const int noutc,
                         const int npeople, const int interaction,
                         const bool iscox)
{
    std::ifstream myfile(fname);
    std::string line;
    std::string tmp;
    noutcomes = noutc;

    int nphenocols = 0;
    int nrpeople = 0;
    if (myfile.is_open())
    {
        /* Read the header and determine the number of columns */
        std::getline(myfile, line);
        std::stringstream line_stream(line);
        // std::cout << line << "\n ";
        while (line_stream >> tmp)
        {
            nphenocols++;
            // std::cout << tmp << " " << nphenocols << " ";
        }

        /* Read the remaining lines */
        while (std::getline(myfile, line))
        {
            int tmplins = 0;
            std::stringstream line_stream(line);
            while (line_stream >> tmp)
                tmplins++;

            if (tmplins != nphenocols)
            {
                std::cerr << "phenofile: number of variables different from "
                          << nphenocols << " in line " << tmplins << endl;
                myfile.close();
                exit(1);
            }
            nrpeople++;
        };
        myfile.close();
    }
    else
    {
        std::cerr << "Unable to open file " << fname << endl;
        exit(1);
    }
    std::cout << "Actual number of people in phenofile = " << nrpeople;

    if (npeople > 0)
    {
        nrpeople = npeople;
        std::cout << "; using only " << nrpeople << " first\n";
    }
    else
    {
        std::cout << "; using all of these\n";
    }

    nids_all = nrpeople;

    /* Determine the number of covariates */
    ncov = nphenocols - 1 - noutcomes;
    model_terms = new std::string[ncov + 2];

    // first pass -- find unmeasured people
    std::ifstream infile(fname);
    if (!infile)
    {
        std::cerr << "phedata: cannot open file " << fname << endl;
    }

    infile >> tmp;
    model = "( ";
    infile >> tmp;
    model = model + tmp;
    for (int i = 1; i < noutcomes; i++)
    {
        infile >> tmp;
        model = model + " , ";
        model = model + tmp;
    }
    n_model_terms = 0;

    model = model + " ) ~ mu";
    model_terms[n_model_terms++] = "mu";

    if (nphenocols > noutcomes + 1) // i.e. we have covariate column(s)
    {
        infile >> tmp;
        model = model + " + " + tmp;
        model_terms[n_model_terms++] = tmp;
        for (int i = (2 + noutcomes); i < nphenocols; i++)
        {
            infile >> tmp;

            model = model + " + ";
            model = model + tmp;
            model_terms[n_model_terms++] = tmp;
        }
    }
    model = model + " + SNP_A1";
    if (interaction != 0)
    {
            model = model + " + "
                + model_terms[interaction]
                + "*SNP_A1";
    }
    model_terms[n_model_terms++] = "SNP_A1";

#if LOGISTIC
    std::cout << "Logistic ";
#elif LINEAR
    std::cout << "Linear ";
#elif COXPH
    std::cout << "Coxph ";
#else
    std::cout << "Unrecognised ";
#endif
    std::cout << "model: " << model << "\n";

    // Filter people with incomplete phenotype (outcome and covariate)
    // data.
    nids = 0;
    for (int i = 0; i < nrpeople; i++)
    {
        /**
         * Add an entry to the allmeasured vector (set to true be
         * default). In the next steps we will determine whether the
         * value of true is correct, or whether this should be set to
         * false because an NA value is encountered.
         */
        allmeasured.push_back(true);

        for (int j = 0; j < nphenocols; j++)
        {
            infile >> tmp;
            if (j > 0 && (tmp[0] == 'N' || tmp[0] == 'n'))
                allmeasured[i] = false;
        }
        if (allmeasured[i])
            nids++;
    }
    infile.close();
    // std::cout << "npeople = " << nids_all
    //           << ", no. all measured = " << nids << "\n";

    // allocate objects
    int ntmpcov = 1;
    if (ncov > 0)
    {
        ntmpcov = ncov;
    }
    idnames = new std::string[nids];

    X.reinit(nids, ntmpcov);
    Y.reinit(nids, noutcomes);
    // second pass -- read the data
    infile.open(fname);
    if (!infile)
    {
        std::cerr << "phedata: cannot open file " << fname << endl;
        exit(1);
    }

    /* Read the header and discard it */
    for (int i = 0; i < nphenocols; i++)
    {
        infile >> tmp;
    }

    /* Fill the X and Y matrices by reading from the pheno file*/
    int m = 0;
    for (int i = 0; i < nrpeople; i++)
        if (allmeasured[i])
        {
            infile >> tmp;
            idnames[m] = tmp;
            for (int j = 0; j < noutcomes; j++)
            {
                infile >> tmp;
                Y.put(std::stod(tmp), m, j);
            }
            for (int j = (1 + noutcomes); j < nphenocols; j++)
            {
                infile >> tmp;
                X.put(std::stod(tmp), m, (j - 1 - noutcomes));
            }
            m++;
        }
        else
        {
            for (int j = 0; j < nphenocols; j++)
                infile >> tmp;
        }
    infile.close();
}


phedata::~phedata()
{
    // delete X;
    // delete Y;
    delete [] model_terms;
    delete [] idnames;
}
