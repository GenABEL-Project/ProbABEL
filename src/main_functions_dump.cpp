/**
 * \file   main_functions_dump.cpp
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Yurii S. Aulchenko (cox, log, lin regressions)
 * \author Maksim V. Struchalin
 *
 * \brief File containing some auxiliary functions that used to be in
 * main.cpp. Having them here helps to keep the overview.
 *
 *  Created on: Nov 27, 2013
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


#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

#include "maskedmatrix.h"
#include "mlinfo.h"
#include "phedata.h"
#include "invsigma.h"
#include "command_line_settings.h"


/**
 * \brief Print progress bar of analysis to stdout.
 *
 * Send a progress update (a percentage) to stdout so that the user
 * has a rough indication of the percentage of SNPs that has already
 * been completed.
 *
 * @param csnp Number of the SNP that is currently being analysed.
 * @param nsnps Total number of SNPs
 */
void update_progress_to_cmd_line(const int csnp, const int nsnps)
{
    std::cout << std::setprecision(2) << std::fixed;

    if (csnp % 1000 == 0)
    {
        if (csnp == 0)
        {
            std::cout << "Analysis: "
                      << std::setw(5)
                      << 100. * static_cast<double>(csnp)
                              / static_cast<double>(nsnps)
                      << "%...";
        }
        else
        {
            std::cout << "\b\b\b\b\b\b\b\b\b"
                      << std::setw(5)
                      << 100. * static_cast<double>(csnp)
                              / static_cast<double>(nsnps)
                      << "%...";
        }
        std::cout.flush();
    }
    std::cout << std::setprecision(6);
}


/**
 * \brief Open an output file for each model when using probability
 * data (npgreds == 2).
 *
 * This function creates the _2df.out.txt etc. files.
 *
 * @param outfile Vector of output streams
 * @param outfilename_str Basename of the outputfiles.
 */
void open_files_for_output(std::vector<std::ofstream*>& outfile,
                           const std::string& outfilename_str)
{
    //create a list of filenames
    const int amount_of_files = 5;
    std::string filenames[amount_of_files] = {
        outfilename_str + "_2df.out.txt",
        outfilename_str + "_add.out.txt",
        outfilename_str + "_domin.out.txt",
        outfilename_str + "_recess.out.txt",
        outfilename_str + "_over_domin.out.txt" };

    for (int i = 0; i < amount_of_files; i++)
    {
        outfile.push_back(new std::ofstream());
        outfile[i]->open((filenames[i]).c_str());
        if (!outfile[i]->is_open())
        {
            std::cerr << "Cannot open file for writing: "
                      << filenames[i]
                      << "\n";
            exit(1);
        }
    }
}


int create_phenotype(phedata& phd, const cmdvars& input_var)
{
    phd.setphedata(input_var.getPhefilename(),
                   input_var.getNoutcomes(),
                   input_var.getNpeople(),
                   input_var.getInteraction(),
                   input_var.isIscox());

    int interaction_cox = input_var.getInteraction();
#if COXPH
    interaction_cox--;
#endif

    if (input_var.getInteraction() < 0 ||
        input_var.getInteraction() > phd.ncov ||
        interaction_cox > phd.ncov)
    {
        std::cerr << "error: Interaction parameter is out of range "
                  << "(interaction="
                  << input_var.getInteraction()
                  << ") \n";
        exit(1);
    }

    return interaction_cox;
}


/**
 * \brief Load the inverse variance-covariance matrix into an InvSigma object.
 *
 * @param input_var Object containing the values of the various
 * command line options.
 * @param phd Object with phenotype data
 * @param invvarmatrix The object of type masked_matrix in which the
 * inverse variance-covariance matrix is returned.
 */
void loadInvSigma(const cmdvars& input_var, const phedata& phd,
                  masked_matrix& invvarmatrix)
{
    std::cout << "You are running mmscore...\n";
    InvSigma inv(input_var.getInverseFilename(), phd);
    // invvarmatrix = inv.get_matrix();
    //double par = 1.; //var(phd.Y)*phd.nids/(phd.nids-phd.ncov-1);
    invvarmatrix.set_matrix(inv.get_matrix());    // = invvarmatrix * par;
    std::cout << " loaded InvSigma...\n" << std::flush;
}


/**
 * \brief Create the first part of the output file header.
 *
 * \param outfile Vector of output file streams. Contains the streams
 * of the output file(s). One file when using dosage data (ngpreds==1)
 * and one for each genetic model in case probabilities are used
 * (ngpreds==2).
 * \param input_var Object containing the values of the various
 * command line options.
 * \param phd Object with phenotype data
 */
void create_start_of_header(std::vector<std::ofstream*>& outfile,
                            const cmdvars& input_var, const phedata& phd)
{
    for (unsigned int i = 0; i < outfile.size(); i++)
    {
        (*outfile[i]) << "name"
                      << input_var.getSep()
                      << "A1"
                      << input_var.getSep()
                      << "A2"
                      << input_var.getSep()
                      << "Freq1"
                      << input_var.getSep()
                      << "MAF"
                      << input_var.getSep()
                      << "Quality"
                      << input_var.getSep()
                      << "Rsq"
                      << input_var.getSep()
                      << "n"
                      << input_var.getSep()
                      << "Mean_predictor_allele";
        if (input_var.getChrom() != "-1")
            (*outfile[i]) << input_var.getSep() << "chrom";
        if (input_var.getMapfilename() != NULL)
            (*outfile[i]) << input_var.getSep() << "position";
        if (input_var.getFlipMAF())
            (*outfile[i]) << input_var.getSep() << "allelesFlipped";
    }

    if (input_var.getAllcov()) //All covariates in output
    {
        for (unsigned int file = 0; file < outfile.size(); file++)
            for (int i = 0; i < phd.n_model_terms - 1; i++)
                *outfile[file] << input_var.getSep()
                               << "beta_"
                               << phd.model_terms[i]
                               << input_var.getSep()
                               << "sebeta_"
                               << phd.model_terms[i];
    }
}


/**
 * \brief Create the rest of header of the output file(s).
 *
 * \sa create_start_of_header for the creation of the first few header
 * columns.
 *
 * \param outfile vector of output file streams. Contains the streams
 * of the output file(s). One file when using dosage data (ngpreds==1)
 * and one for each genetic model in case probabilities are used
 * (ngpreds==2).
 * \param input_var object containing the values of the various
 * command line options.
 * \param phd object with phenotype data
 * \param interaction_cox are we using the Cox model with interaction?
 */
void create_header(std::vector<std::ofstream*>& outfile,
                   const cmdvars& input_var, const phedata& phd,
                   const int& interaction_cox)
{
    create_start_of_header(outfile, input_var, phd);

    if (input_var.getNgpreds() == 1) // dose data: only additive model
    {
        *outfile[0] << input_var.getSep()
                    << "beta_SNP_addA1"
                    << input_var.getSep()
                    << "sebeta_SNP_addA1";

        if (input_var.getInteraction() != 0)
        {
            *outfile[0] << input_var.getSep()
                        << "beta_SNP_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep()
                        << "sebeta_SNP_"
                        << phd.model_terms[interaction_cox];
        }

        if (input_var.getInverseFilename() == NULL)
        {
            //Han Chen
#if !COXPH
            if (input_var.getInteraction() != 0 && !input_var.getAllcov())
            {
                *outfile[0] << input_var.getSep()
                            << "cov_SNP_int_SNP_"
                            << phd.model_terms[interaction_cox];
            }
#endif
        }
        *outfile[0] << input_var.getSep() << "chi2_SNP_add";
        *outfile[0] << "\n";
    } // ngpreds == 1
    else if (input_var.getNgpreds() == 2) // prob data: all models
    {
        *outfile[0] << input_var.getSep()
                    << "beta_SNP_A1A2"
                    << input_var.getSep()
                    << "sebeta_SNP_A1A2"
                    << input_var.getSep()
                    << "beta_SNP_A1A1"
                    << input_var.getSep()
                    << "sebeta_SNP_A1A1";
        *outfile[1] << input_var.getSep()
                    << "beta_SNP_addA1"
                    << input_var.getSep()
                    << "sebeta_SNP_addA1";
        *outfile[2] << input_var.getSep()
                    << "beta_SNP_domA1"
                    << input_var.getSep()
                    << "sebeta_SNP_domA1";
        *outfile[3] << input_var.getSep()
                    << "beta_SNP_recA1"
                    << input_var.getSep()
                    << "sebeta_SNP_recA1";
        *outfile[4] << input_var.getSep()
                    << "beta_SNP_odomA1"
                    << input_var.getSep()
                    << "sebeta_SNP_odomA1";
        if (input_var.getInteraction() != 0)
        {
            //Han Chen
            *outfile[0] << input_var.getSep()
                        << "beta_SNP_A1A2_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep()
                        << "sebeta_SNP_A1A2_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep()
                        << "beta_SNP_A1A1_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep()
                        << "sebeta_SNP_A1A1_"
                        << phd.model_terms[interaction_cox];
#if !COXPH
            if (input_var.getInverseFilename() == NULL &&
                !input_var.getAllcov())
            {
                *outfile[0] << input_var.getSep()
                            << "cov_SNP_A1A2_int_SNP_"
                            << phd.model_terms[interaction_cox]
                            << input_var.getSep()
                            << "cov_SNP_A1A1_int_SNP_"
                            << phd.model_terms[interaction_cox];
            }
#endif
            //Oct 26, 2009
            for (unsigned int file = 1; file < outfile.size(); file++)
            {
                *outfile[file] << input_var.getSep()
                               << "beta_SNP_"
                               << phd.model_terms[interaction_cox]
                               << input_var.getSep()
                               << "sebeta_SNP_"
                               << phd.model_terms[interaction_cox];
                //Han Chen
#if !COXPH
                if (input_var.getInverseFilename() == NULL
                    && !input_var.getAllcov())
                {
                    *outfile[file] << input_var.getSep()
                                   << "cov_SNP_int_SNP_"
                                   << phd.model_terms[interaction_cox];
                }
#endif
                //Oct 26, 2009
            }
        }
        *outfile[0] << input_var.getSep() << "chi2_SNP_2df\n";
        *outfile[1] << input_var.getSep() << "chi2_SNP_add\n";
        *outfile[2] << input_var.getSep() << "chi2_SNP_dom\n";
        *outfile[3] << input_var.getSep() << "chi2_SNP_rec\n";
        *outfile[4] << input_var.getSep() << "chi2_SNP_odom\n";
    } // End: ngpreds == 2
    else
    {
        cerr << "Error: create_header(): ngpreds != 1 or 2.\n";
    }
}


/**
 * \brief Write the information from the mlinfo file to the output
 * file(s).
 *
 * \param outfile Vector of output file(s)
 * \param file index number identifying the file in the vector of files
 * \param mli mlinfo object
 * \param csnp number of the SNP that is currently being analysed
 * \param input_var object containing the information of the options
 * specified on the command line
 * \param gcount The number of non-NaN genotypes
 * \param freq The allele frequency based on the non-NaN genotypes
 */
void write_mlinfo(const std::vector<std::ofstream*>& outfile,
                  const unsigned int file, const mlinfo& mli,
                  const int csnp, const cmdvars& input_var,
                  const int gcount, const double freq)
{
    *outfile[file] << mli.name[csnp]
                   << input_var.getSep()
                   << mli.A1[csnp]
                   << input_var.getSep()
                   << mli.A2[csnp]
                   << input_var.getSep()
                   << mli.Freq1[csnp]
                   << input_var.getSep()
                   << mli.MAF[csnp]
                   << input_var.getSep()
                   << mli.Quality[csnp]
                   << input_var.getSep()
                   << mli.Rsq[csnp]
                   << input_var.getSep()
                   << gcount
                   << input_var.getSep()
                   << freq;
    if (input_var.getChrom() != "-1")
    {
        *outfile[file] << input_var.getSep() << input_var.getChrom();
    }
    if (input_var.getMapfilename() != NULL)
    {
        *outfile[file] << input_var.getSep() << mli.map[csnp];
    }
    if (input_var.getFlipMAF())
    {
        *outfile[file] << input_var.getSep() << mli.allelesFlipped[csnp];
    }
}


/**
 * \brief Get the position within a (row or column) vector (the index)
 * where a \f$ \beta \f$ (or \f$ se_{\beta} \f$) starts.
 *
 * This is basically a matter of counting backwards from the end of
 * the vector/list.
 *
 * @param input_var Object containing the values of the various
 * command line options.
 * @param model Number of the genetic model (additive, etc)
 * @param number_of_rows_or_columns Total number of rows or columns in
 * the vector.
 *
 * @return Start position of beta for this model
 */
int get_start_position(const cmdvars& input_var, const int model,
                       const int number_of_rows_or_columns)
{
    int start_pos;
    if (!input_var.getAllcov() &&
        model == 0 &&
        input_var.getInteraction() == 0)
    {
        if (input_var.getNgpreds() == 2)
        {
            start_pos = number_of_rows_or_columns - 2;
        } else{
            start_pos = number_of_rows_or_columns - 1;
        }
    } else if (!input_var.getAllcov() && model == 0
            && input_var.getInteraction() != 0)
    {
        if (input_var.getNgpreds() == 2)
        {
            start_pos = number_of_rows_or_columns - 4;
        } else{
            start_pos = number_of_rows_or_columns - 2;
        }
    } else if (!input_var.getAllcov() && model != 0
            && input_var.getInteraction() == 0)
    {
        start_pos = number_of_rows_or_columns - 1;
    } else if (!input_var.getAllcov() && model != 0
            && input_var.getInteraction() != 0)
    {
        start_pos = number_of_rows_or_columns - 2;
    } else{
        start_pos = 0;
    }

    return start_pos;
}
