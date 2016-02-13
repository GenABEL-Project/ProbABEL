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


#include <phedata.h>
#include <string>
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

void phedata::set_is_interaction_excluded(const bool int_exl)
{
    is_interaction_excluded = int_exl;
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
    static const unsigned int BFS = 1048576;
    std::ifstream myfile(fname);
    char *line = new char[BFS];
    char *tmp  = new char[BFS];
    std::string interaction_cov_name;
    noutcomes = noutc;

    int nphenocols = 0;
    int nrpeople = 0;
    if (myfile.is_open())
    {
        myfile.getline(line, BFS);
        std::stringstream line_stream(line);
        // std::cout << line << "\n ";
        while (line_stream >> tmp)
        {
            nphenocols++;
            // std::cout << tmp << " " << nphenocols << " ";
        }

        while (myfile.getline(line, BFS))
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

    ncov = nphenocols - 1 - noutcomes;
    nids_all = nrpeople;
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
            if (n_model_terms == interaction && is_interaction_excluded)
               {
                   interaction_cov_name = tmp;
                   n_model_terms++;
                   continue;
               }

            model = model + " + ";
            model = model + tmp;
            model_terms[n_model_terms++] = tmp;
        }
    }
    model = model + " + SNP_A1";
    if (interaction != 0)
    {
        if (iscox)
        {
           if (!is_interaction_excluded)
           {
               model = model + " + "
                   + model_terms[interaction - 1]
                   + "*SNP_A1";
           }
           else
           {
               model = model + " + " + interaction_cov_name + "*SNP_A1";
           }
        }
        else
        {
            if (!is_interaction_excluded)
            {
                model = model + " + " + model_terms[interaction] + "*SNP_A1";
            }
            else
            {
                model = model + " + " + interaction_cov_name + "*SNP_A1";
            }
        }
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
    allmeasured = new unsigned short int[nrpeople];
    nids = 0;
    for (int i = 0; i < nrpeople; i++)
    {
        allmeasured[i] = 1;
        for (int j = 0; j < nphenocols; j++)
        {
            infile >> tmp;
            if (j > 0 && (tmp[0] == 'N' || tmp[0] == 'n'))
                allmeasured[i] = 0;
        }
        if (allmeasured[i] == 1)
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

    for (int i = 0; i < nphenocols; i++)
    {
        infile >> tmp;
    }

    int m = 0;
    for (int i = 0; i < nrpeople; i++)
        if (allmeasured[i] == 1)
        {
            infile >> tmp;
            idnames[m] = tmp;
            for (int j = 0; j < noutcomes; j++)
            {
                infile >> tmp;
                Y.put(atof(tmp), m, j);
            }
            for (int j = (1 + noutcomes); j < nphenocols; j++)
            {
                infile >> tmp;
                X.put(atof(tmp), m, (j - 1 - noutcomes));
            }
            m++;
        }
        else
        {
            for (int j = 0; j < nphenocols; j++)
                infile >> tmp;
        }
    infile.close();

    delete[] line;
    delete[] tmp;
}


phedata::~phedata()
{
    // delete X;
    // delete Y;
    delete [] allmeasured;
    delete [] model_terms;
    delete [] idnames;
}
