/**
 * \file   main.cpp
 * \author Yurii S. Aulchenko (cox, log, lin regressions)
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Maksim V. Struchalin
 *
 * \brief ProbABEL main file
 *
 *
 */
//=============================================================================
//           Filename:  src/main.cpp
//
//        Description:  ProbABEL head file.
//
//             Author:  Yurii S. Aulchenko (cox, log, lin regressions)
//             Modified by: M. Kooyman,
//                          L.C. Karssen,
//                          Maksim V. Struchalin
//
// modified 14-May-2009 by MVS:  interaction with SNP, interaction with SNP
//                               with exclusion of interacted covariates,
//                               mmscore implemented (poor me)
// modified 20-Jul-2009 by YSA: small changes, bug fix in mmscore option
// modified 22-Sep-2009 by YSA: "robust" option added
//
// Modified by Han Chen (hanchen@bu.edu) on Nov 9, 2009
// to extract the covariance between the estimate of beta(SNP) and the estimate
// of beta(interaction) based on src/main.cpp version 0.1-0 as of Oct 19, 2009
//
//  Company:  Department of Epidemiology, ErasmusMC Rotterdam, The Netherlands.
//
//=============================================================================

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


#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

#include <ctime> //needed for timing loading non file vector format

#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#include "maskedmatrix.h"
#include "reg1.h"
#include "command_line_settings.h"
#include "coxph_data.h"
#include "main_functions_dump.h"
#include "mlinfo.h"
#include "invsigma.h"

/**
 * Main routine. The main logic of ProbABEL can be found here
 *
 * \param argc Number of command line arguments
 * \param argv Vector containing the command line arguments
 *
 * \return 0 if all went well. Other integer numbers if an error
 * occurred
 */
int main(int argc, char * argv[])
{
    cmdvars input_var;
    input_var.set_variables(argc, argv);

    input_var.printinfo();

    cout << "Reading info data...\n" << flush;
    mlinfo mli(input_var.getMlinfofilename(),
               input_var.getMapfilename(),
               input_var.getFlipMAF());
    int nsnps = mli.nsnps;
    phedata phd;
    cout << "Reading phenotype data...\n" << flush;
    int interaction_cox = create_phenotype(phd, input_var);

    masked_matrix invvarmatrix;

    if (input_var.getInverseFilename() != NULL)
    {
        loadInvSigma(input_var, phd, invvarmatrix);
    }

    gendata gtd;
    cout << "Reading genotype data... " << flush;
    if (!input_var.getIsFvf())
    {
        // TODO(maartenk): remove timing code
        // make clock to time loading of the non filevector file
        std::clock_t    start;
        start = std::clock();

        // use the non-filevector input format
        gtd.re_gendata(input_var.getGenfilename(), nsnps,
                       input_var.getNgpreds(), phd.nids_all, phd.nids,
                       phd.allmeasured, input_var.getSkipd(), phd.idnames);

        // TODO(maartenk): remove timing code
        double millisec=((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000))/1000;
        cout << "done in "<< millisec<< " seconds.\n" << flush;
    }
    else
    {
        // use the filevector input format (missing second last skipd
        // parameter)
        gtd.re_gendata(input_var.getStrGenfilename(), nsnps,
                       input_var.getNgpreds(), phd.nids_all, phd.nids,
                       phd.allmeasured, phd.idnames);
        cout << "done.\n" << flush;
    }


    // Set the number and names of the genetic models that will be
    // analysed.
    // The total number of genetic models
    int maxmod = 0;
    // The human readable names of the genetic models
    std::vector<std::string> modelNames;

    if (input_var.getNgpreds() == 1) {
        // For dosage data we can only run the additive model.
        maxmod = 1;
        modelNames.push_back("additive");
    } else if (input_var.getNgpreds() == 2) {
        maxmod = 5;             // Only with probability data can we
                                // run all of them.
        modelNames.push_back("2df");
        modelNames.push_back("additive");
        modelNames.push_back("dominant");
        modelNames.push_back("recessive");
        modelNames.push_back("overdominant");
    }


    // estimate null model
#if COXPH
    coxph_data nrgd = coxph_data(phd, gtd, -1);
#else
    regdata nrgd = regdata(phd, gtd, -1);
#endif

    std::cout << " loaded null data..." << std::flush;
#if LOGISTIC
    logistic_reg nrd = logistic_reg(nrgd);

    nrd.estimate(0, 0,
                 input_var.getInteraction(),
                 input_var.getNgpreds(),
                 invvarmatrix,
                 input_var.getRobust(),
                 1);
#elif LINEAR

    linear_reg nrd = linear_reg(nrgd);
#if DEBUG
    std::cout << "[DEBUG] linear_reg nrd = linear_reg(nrgd); DONE.";
#endif
    nrd.estimate(0, 0, input_var.getInteraction(),
                 input_var.getNgpreds(), invvarmatrix,
                 input_var.getRobust(), 1);
#elif COXPH
    coxph_reg nrd = coxph_reg(nrgd);
    nrd.estimate(nrgd, 0, modelNames,
                 input_var.getInteraction(), input_var.getNgpreds(),
                 true, 1, mli, 0);
#endif
    double null_loglik = nrd.loglik;

    std::cout << " estimated null model...";
    // end null
#if COXPH
    coxph_data rgd(phd, gtd, 0);
#else
    regdata rgd(phd, gtd, 0);
#endif
    std::cout << " formed regression object...\n";


    // Open a vector of files that will be used for output. Depending
    // on the number of genomic predictors we either open 5 files (one
    // for each model if we have prob data) or one (if we have dosage
    // data).
    std::string outfilename_str(input_var.getOutfilename());
    std::vector<std::ofstream*> outfile;

    // Prob data: All models output. One file per model
    if (input_var.getNgpreds() == 2)
    {
        open_files_for_output(outfile, outfilename_str);
        if (input_var.getNohead() != 1)
        {
            create_header(outfile, input_var, phd, interaction_cox);
        }
    }
    else  // Dosage data: Only additive model => only one output file
    {
        outfile.push_back(
            new std::ofstream((outfilename_str + "_add.out.txt").c_str()));

        if (!outfile[0]->is_open())
        {
            std::cerr << "Cannot open file for writing: "
                      << outfilename_str
                      << "\n";
            exit(1);
        }
        if (input_var.getNohead() != 1)
        {
            create_header(outfile, input_var, phd, interaction_cox);
        }
    }  // END else: we have dosage data => only one file

    int start_pos, end_pos;

    std::vector<std::ostringstream *> beta_sebeta;
    // Han Chen
    std::vector<std::ostringstream *> covvalue;
    // Oct 26, 2009
    std::vector<std::ostringstream *> chi2;

    // Create string streams for betas, SEs, etc. These are used to
    // later store the various output values that will be written to
    // files.
    for (int i = 0; i < maxmod; i++)
    {
        beta_sebeta.push_back(new std::ostringstream());
        beta_sebeta[i]->precision(6);
        // *beta_sebeta[i] << scientific;
        // Han Chen
        covvalue.push_back(new std::ostringstream());
        covvalue[i]->precision(6);
        // *covvalue[i] << scientific;
        // Oct 26, 2009
        chi2.push_back(new std::ostringstream());
        chi2[i]->precision(6);
        // *chi2[i] << scientific;
    }


    // Here we start the analysis for each SNP.
    for (int csnp = 0; csnp < nsnps; csnp++)
    {
        rgd.update_snp(&gtd, csnp, mli, input_var.getFlipMAF());


        int poly = 1;
        if (fabs(rgd.freq) < 1.e-16 || fabs(1. - rgd.freq) < 1.e-16)
        {
            poly = 0;
        }

        if (fabs(mli.Rsq[csnp]) < 1.e-16)
        {
            poly = 0;
        }

        // Write mlinfo information to the output file(s)
        // Prob data: All models output. One file per model
        if (input_var.getNgpreds() == 2)
        {
            for (unsigned int file = 0; file < outfile.size(); file++)
            {
                write_mlinfo(outfile, file, mli, csnp, input_var,
                             rgd.gcount, rgd.freq);
            }
        } else{
            // Dosage data: only additive model
            int file = 0;
            write_mlinfo(outfile, file, mli, csnp, input_var,
                         rgd.gcount, rgd.freq);
        }

        // Run regression for each model for the current SNP
        for (int model = 0; model < maxmod; model++)
        {
            if (poly) // Allele freq is not too rare
            {
#if LOGISTIC
                logistic_reg rd(rgd);
#elif LINEAR
                linear_reg rd(rgd);
#elif COXPH
                coxph_reg rd(rgd);
#endif
#if !COXPH
                if (input_var.getScore())
                {
                    rd.score(nrd.residuals, model,
                             input_var.getInteraction(),
                             input_var.getNgpreds(),
                             invvarmatrix);
                }
                else
                {
                    rd.estimate(0, model,
                                input_var.getInteraction(),
                                input_var.getNgpreds(),
                                invvarmatrix,
                                input_var.getRobust());
                }
#else
                rd.estimate(rgd, model,
                            modelNames,
                            input_var.getInteraction(),
                            input_var.getNgpreds(), true, 0, mli, csnp);
#endif

                int number_of_rows_or_columns = rd.beta.nrow;
                start_pos = get_start_position(input_var, model,
                                               number_of_rows_or_columns);

                // The regression coefficients for the SNPs are in the
                // last rows of beta[] and sebeta[].
                for (int pos = start_pos; pos < rd.beta.nrow; pos++)
                {
                    *beta_sebeta[model] << input_var.getSep()
                                        << rd.beta[pos]
                                        << input_var.getSep()
                                        << rd.sebeta[pos];
                    // Han Chen
#if !COXPH
                    if (input_var.getInverseFilename() == NULL
                            && !input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        if (pos > start_pos)
                        {
                            if (model == 0)
                            {
                                if (input_var.getNgpreds() == 2)
                                {
                                    if (pos > start_pos + 2)
                                    {
                                        *covvalue[model]
                                            << rd.covariance[pos - 3]
                                            << input_var.getSep()
                                            << rd.covariance[pos - 2];
                                    }
                                }  // END ngpreds=2
                                else
                                {
                                    *covvalue[model] << rd.covariance[pos - 1];
                                }
                            }  // END model == 0
                            else
                            {
                                *covvalue[model] << rd.covariance[pos - 1];
                            }  // END model != 0
                        }  // END if pos > start_pos
                    }
#endif
                    // Oct 26, 2009
                }  // END for(pos = start_pos; pos < rd.beta.nrow; pos++)


                // calculate chi^2
                // ________________________________
                // cout <<  rd.loglik<<" "<<input_var.getNgpreds() << "\n";

                if (input_var.getInverseFilename() == NULL)
                { // Only if we don't have an inv.sigma file can we use LRT
                    if (input_var.getScore() == 0)
                    {
                        double loglik = rd.loglik;
                        if (rgd.gcount != gtd.nids)
                        {
                            // If SNP data is missing we didn't
                            // correctly compute the null likelihood

                            // Recalculate null likelihood by
                            // stripping the SNP data column(s) from
                            // the X matrix in the regression object
                            // and run the null model estimation again
                            // for this SNP.
#if !COXPH
                            regdata new_rgd = rgd;
#else
                            coxph_data new_rgd = rgd;
#endif

                            new_rgd.remove_snp_from_X();

#ifdef LINEAR
                            linear_reg new_null_rd(new_rgd);
#elif LOGISTIC
                            logistic_reg new_null_rd(new_rgd);
#endif
#if !COXPH
                            new_null_rd.estimate(0,
                                                 model,
                                                 input_var.getInteraction(),
                                                 input_var.getNgpreds(),
                                                 invvarmatrix,
                                                 input_var.getRobust(), 1);
#else
                            coxph_reg new_null_rd(new_rgd);
                            new_null_rd.estimate(new_rgd,
                                                 model, modelNames,
                                                 input_var.getInteraction(),
                                                 input_var.getNgpreds(),
                                                 true, 1, mli, csnp);
#endif
                            *chi2[model] << 2. * (loglik - new_null_rd.loglik);
                        }
                        else
                        {
                            // No missing SNP data, we can compute the LRT
                            *chi2[model] << 2. * (loglik - null_loglik);
                        }
                    } else{
                        // We want score test output
                        *chi2[model] << rd.chi2_score;
                    }
                }  // END if( inv.sigma == NULL )
                else if (input_var.getInverseFilename() != NULL)
                {
                    // We can't use the LRT here, because mmscore is a
                    // REML method. Therefore go for the Wald test
                    if (input_var.getNgpreds() == 2 && model == 0)
                    {
                        /* For the 2df model we can't simply use the
                         * Wald statistic. This can be fixed using the
                         * equation just below Eq.(4) in the ProbABEL
                         * paper. TODO LCK
                         */
                        *chi2[model] << NAN;
                    }
                    else
                    {
                        double Z = rd.beta[start_pos] / rd.sebeta[start_pos];
                        *chi2[model] << Z * Z;
                    }
                }
            }  // END first part of if(poly); allele not too rare
            else
            {   // SNP is rare: beta, sebeta, chi2 = NaN

                // Find out how many columns of nan need to be
                // printed.
                end_pos = 1;    // One entry for the SNP term

                if (input_var.getAllcov())
                {
                    end_pos += phd.n_model_terms - 1;
                }

                if(input_var.getNgpreds() == 2 && model == 0)
                {
                    // The 2df model needs an extra entry for the
                    // second genotype coefficient.
                    end_pos++;
                }

                if(input_var.getInteraction() != 0)
                {
                    // There's an interaction term, add an extra entry
                    // for its coefficient.
                    end_pos++;
                    if(input_var.getNgpreds() == 2 && model == 0)
                    {
                        // For interaction + 2df model another entry
                        // needs to be added, because each of the
                        // genotype coefficients for the interaction.
                        end_pos++;
                    }
                }

                for (int pos = 0; pos < end_pos; pos++)
                {
                    *beta_sebeta[model] << input_var.getSep()
                            << NAN
                            << input_var.getSep()
                            << NAN;
                }

                if (input_var.getNgpreds() == 2)
                {
                    // Han Chen
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        if (model == 0)
                        {
                            *covvalue[model] << NAN
                                             << input_var.getSep()
                                             << NAN;
                        } else{
                            *covvalue[model] << NAN;
                        }
                    }
#endif
                    // Oct 26, 2009
                    *chi2[model] << NAN;
                } else{
                    // ngpreds==1 (and SNP is rare)
                    if (input_var.getInverseFilename() == NULL)
                    {
                        //                     Han Chen
#if !COXPH
                        if (!input_var.getAllcov()
                                && input_var.getInteraction() != 0)
                        {
                            *covvalue[model] << NAN;
                        }
#endif
                        // Oct 26, 2009
                    }  // END if getInverseFilename == NULL
                    *chi2[model] << NAN;
                }  // END ngpreds == 1 (and SNP is rare)
            }  // END else: SNP is rare
        }  // END of model cycle


        // Start writing beta's, se_beta's etc. to file
        if (input_var.getNgpreds() == 2)
        {
            for (int model = 0; model < maxmod; model++)
            {
                *outfile[model] << beta_sebeta[model]->str()
                                << input_var.getSep();
#if !COXPH
                if (!input_var.getAllcov() && input_var.getInteraction() != 0)
                {
                    *outfile[model] << covvalue[model]->str()
                                    << input_var.getSep();
                }
#endif
                *outfile[model] << chi2[model]->str()
                                << "\n";
            }  // END for loop over all models
        }
        else  // Dose data: only additive model. Only one output file
        {
            *outfile[0] << beta_sebeta[0]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[0] << covvalue[0]->str() << input_var.getSep();
            }
#endif
            *outfile[0] << chi2[0]->str() << "\n";
        }  // End ngpreds == 1 when writing output files


        // Clean chi2 and other streams
        for (int model = 0; model < maxmod; model++)
        {
            beta_sebeta[model]->str("");
            // Han Chen
            covvalue[model]->str("");
            // Oct 26, 2009
            chi2[model]->str("");
        }

        update_progress_to_cmd_line(csnp, nsnps);
    }  // END for loop over all SNPs


    // We're almost done. All computations have finished, time to
    // clean up.

    std::cout << setprecision(2) << fixed;
    std::cout << "\b\b\b\b\b\b\b\b\b" << 100.;
    std::cout << "%... done\n";

    // Close output files
    for (unsigned int i = 0; i < outfile.size(); i++)
    {
        outfile[i]->close();
        delete outfile[i];
    }

    // delete gtd;

    // Clean up a couple of vectors
    std::vector<std::ostringstream *>::iterator it = beta_sebeta.begin();
    while (it != beta_sebeta.end())
    {
        delete *it;
        ++it;
    }
    it = covvalue.begin();
    while (it != covvalue.end())
    {
        delete *it;
        ++it;
    }
    it = chi2.begin();
    while (it != chi2.end())
    {
        delete *it;
        ++it;
    }

    return (0);
}
