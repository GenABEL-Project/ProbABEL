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
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif
#include "maskedmatrix.h"
#include "data.h"
#include "reg1.h"
#include "command_line_settings.h"

#include "coxph_data.h"
//#include "coxph_reg.cpp"

#define MAXITER 10

#define EPS 1.e-8
#define CHOLTOL 1.5e-12

void update_progress_to_cmd_line(int csnp, int nsnps)
{
    std::cout << setprecision(2) << fixed;

    if (csnp % 1000 == 0)
    {
        if (csnp == 0)
        {
            std::cout << "Analysis: " << setw(5)
                    << 100. * static_cast<double>(csnp)
                            / static_cast<double>(nsnps) << "%...";
        } else
        {
            `
            std::cout << "\b\b\b\b\b\b\b\b\b" << setw(5)
                    << 100. * static_cast<double>(csnp)
                            / static_cast<double>(nsnps) << "%...";
        }
        std::cout.flush();
    }
}

void open_files_for_output(std::vector<std::ofstream*>& outfile,
        std::string& outfilename_str)
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
            std::cerr << "Cannot open file for writing: " << filenames[i]
                    << "\n";
            exit(1);
        }
    }
}

int create_phenotype(phedata& phd, cmdvars& input_var)
{
    phd.set_is_interaction_excluded(input_var.isIsInteractionExcluded());
    phd.setphedata(input_var.getPhefilename(), input_var.getNoutcomes(),
            input_var.getNpeople(), input_var.getInteraction(),
            input_var.isIscox());

    int interaction_cox = input_var.getInteraction();
#if COXPH
    interaction_cox--;
#endif

    if (input_var.getInteraction() < 0 || input_var.getInteraction() > phd.ncov
            || interaction_cox > phd.ncov)
    {
        std::cerr << "error: Interaction parameter is out of range "
                << "(interaction=" << input_var.getInteraction() << ") \n";
        exit(1);
    }

    return interaction_cox;
}

void loadInvSigma(cmdvars& input_var, phedata& phd, masked_matrix& invvarmatrix)
{
    std::cout << "you are running mmscore...\n";
    InvSigma inv(input_var.getInverseFilename(), &phd);
    // invvarmatrix = inv.get_matrix();
    //double par = 1.; //var(phd.Y)*phd.nids/(phd.nids-phd.ncov-1);
    invvarmatrix.set_matrix(inv.get_matrix());    // = invvarmatrix * par;
    std::cout << " loaded InvSigma ..." << std::flush;
}

void create_start_of_header(std::vector<std::ofstream*>& outfile,
        cmdvars& input_var, phedata& phd)
{
    for (unsigned int i = 0; i < outfile.size(); i++)
    {
        (*outfile[i]) << "name" << input_var.getSep() << "A1"
                << input_var.getSep() << "A2" << input_var.getSep() << "Freq1"
                << input_var.getSep() << "MAF" << input_var.getSep()
                << "Quality" << input_var.getSep() << "Rsq"
                << input_var.getSep() << "n" << input_var.getSep()
                << "Mean_predictor_allele";
        if (input_var.getChrom() != "-1")
            (*outfile[i]) << input_var.getSep() << "chrom";
        if (input_var.getMapfilename() != NULL)
            (*outfile[i]) << input_var.getSep() << "position";
    }

    if (input_var.getAllcov()) //All covariates in output
    {
        for (unsigned int file = 0; file < outfile.size(); file++)
            for (int i = 0; i < phd.n_model_terms - 1; i++)
                *outfile[file] << input_var.getSep() << "beta_"
                        << phd.model_terms[i] << input_var.getSep() << "sebeta_"
                        << phd.model_terms[i];
    }
}

void create_header_1(std::vector<std::ofstream*>& outfile, cmdvars& input_var,
        phedata& phd, int& interaction_cox)
{
    create_start_of_header(outfile, input_var, phd);

    *outfile[0] << input_var.getSep() << "beta_SNP_A1A2" << input_var.getSep()
            << "beta_SNP_A1A1" << input_var.getSep() << "sebeta_SNP_A1A2"
            << input_var.getSep() << "sebeta_SNP_A1A1";

    *outfile[1] << input_var.getSep() << "beta_SNP_addA1" << input_var.getSep()
            << "sebeta_SNP_addA1";
    *outfile[2] << input_var.getSep() << "beta_SNP_domA1" << input_var.getSep()
            << "sebeta_SNP_domA1";
    *outfile[3] << input_var.getSep() << "beta_SNP_recA1" << input_var.getSep()
            << "sebeta_SNP_recA1";
    *outfile[4] << input_var.getSep() << "beta_SNP_odom" << input_var.getSep()
            << "sebeta_SNP_odom";
    if (input_var.getInteraction() != 0)
    {
        //Han Chen
        *outfile[0] << input_var.getSep() << "beta_SNP_A1A2_"
                << phd.model_terms[interaction_cox] << input_var.getSep()
                << "sebeta_SNP_A1A2_" << phd.model_terms[interaction_cox]
                << input_var.getSep() << "beta_SNP_A1A1_"
                << phd.model_terms[interaction_cox] << input_var.getSep()
                << "sebeta_SNP_A1A1_" << phd.model_terms[interaction_cox];
#if !COXPH
        if (input_var.getInverseFilename() == NULL && !input_var.getAllcov())
        {
            *outfile[0] << input_var.getSep() << "cov_SNP_A1A2_int_SNP_"
                    << phd.model_terms[interaction_cox] << input_var.getSep()
                    << "cov_SNP_A1A1_int_SNP_"
                    << phd.model_terms[interaction_cox];
        }
#endif
        //Oct 26, 2009
        for (unsigned int file = 1; file < outfile.size(); file++)
        {
            *outfile[file] << input_var.getSep() << "beta_SNP_"
                    << phd.model_terms[interaction_cox] << input_var.getSep()
                    << "sebeta_SNP_" << phd.model_terms[interaction_cox];
            //Han Chen
#if !COXPH
            if (input_var.getInverseFilename() == NULL
                    && !input_var.getAllcov())
            {
                *outfile[file] << input_var.getSep() << "cov_SNP_int_SNP_"
                        << phd.model_terms[interaction_cox];
            }
#endif
            //Oct 26, 2009
        }
    }
    *outfile[0] << input_var.getSep() << "loglik\n"; //"chi2_SNP_2df\n";
    *outfile[1] << input_var.getSep() << "loglik\n"; //"chi2_SNP_A1\n";
    *outfile[2] << input_var.getSep() << "loglik\n"; //"chi2_SNP_domA1\n";
    *outfile[3] << input_var.getSep() << "loglik\n"; //"chi2_SNP_recA1\n";
    *outfile[4] << input_var.getSep() << "loglik\n"; //"chi2_SNP_odom\n";
}

void create_header2(std::vector<std::ofstream*>& outfile, cmdvars& input_var,
        phedata& phd, int interaction_cox)
{
    create_start_of_header(outfile, input_var, phd);
    *outfile[0] << input_var.getSep() << "beta_SNP_add" << input_var.getSep()
            << "sebeta_SNP_add";

    if (input_var.getInteraction() != 0)
    {
        *outfile[0] << input_var.getSep() << "beta_SNP_"
                << phd.model_terms[interaction_cox] << input_var.getSep()
                << "sebeta_SNP_" << phd.model_terms[interaction_cox];
    }

    if (input_var.getInverseFilename() == NULL)
    {
        //Han Chen
#if !COXPH
        if (input_var.getInteraction() != 0 && !input_var.getAllcov())
        {
            *outfile[0] << input_var.getSep() << "cov_SNP_int_SNP_"
                    << phd.model_terms[interaction_cox];
        }
#endif
        *outfile[0] << input_var.getSep() << "loglik"; //"chi2_SNP";
    }
    //Oct 26, 2009
    *outfile[0] << "\n";
}

void write_mlinfo(const std::vector<std::ofstream*>& outfile, unsigned int file,
        const mlinfo& mli, int csnp, const cmdvars& input_var, int gcount,
        double freq)
{
    *outfile[file] << mli.name[csnp] << input_var.getSep() << mli.A1[csnp]
            << input_var.getSep() << mli.A2[csnp] << input_var.getSep()
            << mli.Freq1[csnp] << input_var.getSep() << mli.MAF[csnp]
            << input_var.getSep() << mli.Quality[csnp] << input_var.getSep()
            << mli.Rsq[csnp] << input_var.getSep() << gcount
            << input_var.getSep() << freq;
    if (input_var.getChrom() != "-1")
    {
        *outfile[file] << input_var.getSep() << input_var.getChrom();
    }
    if (input_var.getMapfilename() != NULL)
    {
        *outfile[file] << input_var.getSep() << mli.map[csnp];
    }
}

int main(int argc, char * argv[])
{
    cmdvars input_var;
    input_var.set_variables(argc, argv);

    input_var.printinfo();
    //	if (allcov && ngpreds>1)
    //	{
    //      cout << "\n\n"
    //           << "WARNING: --allcov allowed only for 1 predictor (MLDOSE)\n";
    //		allcov = 0;
    //	}
    mlinfo mli(input_var.getMlinfofilename(), input_var.getMapfilename());
    int nsnps = mli.nsnps;
    phedata phd;
    int interaction_cox = create_phenotype(phd, input_var);

    //interaction--;
    //	if (input_var.getInverseFilename()!= NULL && phd.ncov > 1)
    //     {
    //         std::cerr << "Error: In mmscore you can not use any covariates."
    //                   << " You phenotype file must conatin id column and "
    //                   << "trait (residuals) only\n";
    //         exit(1);
    //      }
    //	if (input_var.getInverseFilename()!= NULL &&
    //      (allcov == 1 || score == 1
    //                   || input_var.getInteraction()!= 0
    //                   || ngpreds==2))
    //      {
    //          std::cerr << "Error: In mmscore you can use additive model "
    //                    << "without any inetractions only\n";
    //          exit(1);
    //      }
    masked_matrix invvarmatrix;

    /*
     * now should be possible... delete this part later when everything works
     #if LOGISTIC
     if (input_var.getInverseFilename()!= NULL)
     {
     std::cerr << "ERROR: mmscore is forbidden for logistic regression\n";
     exit(1);
     }
     #endif
     */

    std::cout << "Reading data ..." << std::flush;
    if (input_var.getInverseFilename() != NULL)
    {
        loadInvSigma(input_var, phd, invvarmatrix);
    }

    gendata gtd;
    if (!input_var.getIsFvf())
        // use the non-filevector input format
        gtd.re_gendata(input_var.getGenfilename(), nsnps,
                input_var.getNgpreds(), phd.nids_all, phd.nids, phd.allmeasured,
                input_var.getSkipd(), phd.idnames);
    else
        // use the filevector input format (missing second last skipd
        // parameter)
        gtd.re_gendata(input_var.getStrGenfilename(), nsnps,
                input_var.getNgpreds(), phd.nids_all, phd.nids, phd.allmeasured,
                phd.idnames);

    std::cout << " loaded genotypic data ..." << std::flush;
    /**
     if (input_var.getIsFvf())
     gendata gtd(str_genfilename, nsnps, input_var.getNgpreds(),
     phd.nids_all, phd.allmeasured, phd.idnames);
     else
     gendata gtd(input_var.getGenfilename(), nsnps,
     input_var.getNgpreds(), phd.nids_all, phd.nids,
     phd.allmeasured, skipd, phd.idnames);
     **/

    // estimate null model
#if COXPH
    coxph_data nrgd = coxph_data(phd, gtd, -1);
#else
    regdata nrgd = regdata(phd, gtd, -1, input_var.isIsInteractionExcluded());
#endif

    std::cout << " loaded null data ..." << std::flush;
#if LOGISTIC
    logistic_reg nrd = logistic_reg(nrgd);
    nrd.estimate(nrgd, 0, MAXITER, EPS, CHOLTOL, 0,
            input_var.getInteraction(), input_var.getNgpreds(),
            invvarmatrix, input_var.getRobust(), 1);
#elif LINEAR

    linear_reg nrd = linear_reg(nrgd);
#if DEBUG
    std::cout << "[DEBUG] linear_reg nrd = linear_reg(nrgd); DONE.";
#endif
    nrd.estimate(nrgd, 0, CHOLTOL, 0, input_var.getInteraction(),
            input_var.getNgpreds(), invvarmatrix, input_var.getRobust(), 1);
#elif COXPH
    coxph_reg nrd(nrgd);

    nrd.estimate(nrgd, 0, MAXITER, EPS, CHOLTOL, 0,
            input_var.getInteraction(), input_var.getNgpreds(), 1);
#endif

    std::cout << " estimated null model ...";
    // end null
#if COXPH
    coxph_data rgd(phd, gtd, 0);
#else
    regdata rgd(phd, gtd, 0, input_var.isIsInteractionExcluded());
#endif

    std::cout << " formed regression object ...";
    std::cout << " done\n" << std::flush;

    //________________________________________________________________
    //Maksim, 9 Jan, 2009
    std::string outfilename_str(input_var.getOutfilename());
    std::vector<std::ofstream*> outfile;

    //All models output.One file per each model
    if (input_var.getNgpreds() == 2)
    {
        open_files_for_output(outfile, outfilename_str);
        if (input_var.getNohead() != 1)
        {
            create_header_1(outfile, input_var, phd, interaction_cox);
        }
    } else //Only additive model. Only one output file
    {
        outfile.push_back(
                new std::ofstream((outfilename_str + "_add.out.txt").c_str()));

        if (!outfile[0]->is_open())
        {
            std::cerr << "Cannot open file for writing: " << outfilename_str
                    << "\n";
            exit(1);
        }
        if (input_var.getNohead() != 1)
        {
            create_header2(outfile, input_var, phd, interaction_cox);
        }
    }

    //________________________________________________________________

    /*
     if (input_var.getAllcov())
     {
     if (score)
     {
     outfile << input_var.getSep() << "beta_mu"; // << input_var.getSep() << "beta_SNP_A1";
     outfile << input_var.getSep() << "sebeta_mu"; // << input_var.getSep() << "sebeta_SNP_A1";
     }
     else
     {
     for (int i =0; i<phd.n_model_terms-1;i++)
     outfile << input_var.getSep() << "beta_" << phd.model_terms[i] << input_var.getSep() << "sebeta_" << phd.model_terms[i];
     }
     if(interactio != 0) outfile << input_var.getSep() << "beta_SNP_" << phd.model_terms[interaction];
     }
     if (input_var.getNgpreds()==2)
     {
     outfile << input_var.getSep() << "beta_SNP_A1A2"
     << input_var.getSep() << "beta_SNP_A1A1"
     << input_var.getSep() << "sebeta_SNP_A1A2"
     << input_var.getSep() << "sebeta_SNP_a1A1"
     << input_var.getSep() << "chi2_SNP_2df"
     << input_var.getSep() << "beta_SNP_addA1"
     << input_var.getSep() << "sebeta_SNP_addA1"
     << input_var.getSep() << "chi2_SNP_addA1"
     << input_var.getSep() << "beta_SNP_domA1"
     << input_var.getSep() << "sebeta_SNP_domA1"
     << input_var.getSep() << "chi2_SNP_domA1"
     << input_var.getSep() << "beta_SNP_recA1"
     << input_var.getSep() << "sebeta_SNP_recA1"
     << input_var.getSep() << "chi2_SNP_recA1"
     << input_var.getSep() << "beta_SNP_odom"
     << input_var.getSep() << "sebeta_SNP_odom"
     << input_var.getSep() << "chi2_SNP_odom\n";
     }
     else
     {
     outfile << input_var.getSep() << "beta_SNP_add"
     << input_var.getSep() << "sebeta_SNP_add"
     << input_var.getSep() << "chi2_SNP_add\n";
     }

     */
    //	exit(1);
    //________________________________________________________________
    //Maksim, 9 Jan, 2009
    int maxmod = 5;
    int start_pos, end_pos;

    std::vector<std::ostringstream *> beta_sebeta;
    //Han Chen
    std::vector<std::ostringstream *> covvalue;
    //Oct 26, 2009
    std::vector<std::ostringstream *> chi2;

    for (int i = 0; i < maxmod; i++)
    {
        beta_sebeta.push_back(new std::ostringstream());
        //Han Chen
        covvalue.push_back(new std::ostringstream());
        //Oct 26, 2009
        chi2.push_back(new std::ostringstream());
    }

    for (int csnp = 0; csnp < nsnps; csnp++)
    {
        rgd.update_snp(gtd, csnp);
        double freq = 0.;
        int gcount = 0;
        float snpdata1[gtd.nids];
        float snpdata2[gtd.nids];
        if (input_var.getNgpreds() == 2)
        {
            //freq = ((gtd.G).column_mean(csnp*2)*2. +
            //        (gtd.G).column_mean(csnp*2+1))/2.;
            gtd.get_var(csnp * 2, snpdata1);
            gtd.get_var(csnp * 2 + 1, snpdata2);
            for (unsigned int ii = 0; ii < gtd.nids; ii++)
                if (!isnan(snpdata1[ii]) && !isnan(snpdata2[ii]))
                {
                    gcount++;
                    freq += snpdata1[ii] + snpdata2[ii] * 0.5;
                }
        } else
        {
            // freq = (gtd.G).column_mean(csnp)/2.;
            gtd.get_var(csnp, snpdata1);
            for (unsigned int ii = 0; ii < gtd.nids; ii++)
                if (!isnan(snpdata1[ii]))
                {
                    gcount++;
                    freq += snpdata1[ii] * 0.5;
                }
        }
        freq /= static_cast<double>(gcount);
        int poly = 1;
        if (fabs(freq) < 1.e-16 || fabs(1. - freq) < 1.e-16)
        {
            poly = 0;
        }
        if (fabs(mli.Rsq[csnp]) < 1.e-16)
        {
            poly = 0;
        }
        //
        //All models output. One file per each model
        //
        if (input_var.getNgpreds() == 2)
        {
            //Write mlinfo to output:
            for (unsigned int file = 0; file < outfile.size(); file++)
            {
                write_mlinfo(outfile, file, mli, csnp, input_var, gcount, freq);
            }

            for (int model = 0; model < maxmod; model++)
            {
                if (poly) //allel freq is not to rare
                {

#if LOGISTIC
                    logistic_reg rd(rgd);
                    if (input_var.getScore())
                    rd.score(nrd.residuals, rgd, 0, CHOLTOL, model,
                            input_var.getInteraction(),
                            input_var.getNgpreds(),
                            invvarmatrix);
                    else
                    rd.estimate(rgd, 0, MAXITER, EPS, CHOLTOL, model,
                            input_var.getInteraction(),
                            input_var.getNgpreds(),
                            invvarmatrix,
                            input_var.getRobust());
#elif LINEAR
                    linear_reg rd(rgd);
                    if (input_var.getScore())
                        rd.score(nrd.residuals, rgd, 0, CHOLTOL, model,
                                input_var.getInteraction(),
                                input_var.getNgpreds(), invvarmatrix);
                    else
                    {
                        //	rd.mmscore(rgd,0,CHOLTOL,model,input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix);
                        rd.estimate(rgd, 0, CHOLTOL, model,
                                input_var.getInteraction(),
                                input_var.getNgpreds(), invvarmatrix,
                                input_var.getRobust());
                    }
#elif COXPH
                    coxph_reg rd(rgd);
                    rd.estimate(rgd, 0, MAXITER, EPS, CHOLTOL, model,
                            input_var.getInteraction(), true,
                            input_var.getNgpreds());
#endif

                    if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() == 0)
                    {
                        if (input_var.getNgpreds() == 2)
                        {
                            start_pos = rd.beta.nrow - 2;
                        } else
                        {
                            start_pos = rd.beta.nrow - 1;
                        }
                    } else if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() != 0)
                    {
                        if (input_var.getNgpreds() == 2)
                        {
                            start_pos = rd.beta.nrow - 4;
                        } else
                        {
                            start_pos = rd.beta.nrow - 2;
                        }
                    } else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() == 0)
                    {
                        start_pos = rd.beta.nrow - 1;
                    } else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() != 0)
                    {
                        start_pos = rd.beta.nrow - 2;
                    } else
                    {
                        start_pos = 0;
                    }

                    for (int pos = start_pos; pos < rd.beta.nrow; pos++)
                    {
                        *beta_sebeta[model] << input_var.getSep()
                                << rd.beta[pos] << input_var.getSep()
                                << rd.sebeta[pos];
                        //Han Chen
#if !COXPH
                        if (input_var.getInverseFilename() == NULL
                                && !input_var.getAllcov()
                                && input_var.getInteraction() != 0)
                        {
                            if (pos > start_pos)
                            {
                                if (model == 0)
                                {
                                    if (pos > start_pos + 2
                                            || input_var.getNgpreds() != 2)
                                    {
                                        *covvalue[model]
                                                << rd.covariance[pos - 3]
                                                << input_var.getSep()
                                                << rd.covariance[pos - 2];
                                    }
                                } else
                                {
                                    *covvalue[model] << rd.covariance[pos - 1];
                                }
                            }
                        }
#endif
                        //Oct 26, 2009
                    }

                    //calculate chi2
                    //________________________________
                    if (input_var.getInverseFilename() == NULL
                            && input_var.getNgpreds() != 2)
                    {

                        if (input_var.getScore() == 0)
                        {
                            //*chi2[model] << 2.*(rd.loglik-null_loglik);
                            *chi2[model] << rd.loglik;
                        } else
                        {
                            //*chi2[model] << rd.chi2_score;
                            *chi2[model] << "nan";
                        }
                    }
                    //________________________________
                } else //beta, sebeta = nan
                {
                    if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() == 0)
                    {
                        if (input_var.getNgpreds() == 2)
                        {
                            start_pos = rgd.X.ncol - 2;
                        } else
                        {
                            start_pos = rgd.X.ncol - 1;
                        }
                    } else if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() != 0)
                    {
                        if (input_var.getNgpreds() == 2)
                        {
                            start_pos = rgd.X.ncol - 4;
                        } else
                        {
                            start_pos = rgd.X.ncol - 2;
                        }
                    } else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() == 0)
                    {
                        start_pos = rgd.X.ncol - 1;
                    } else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() != 0)
                    {
                        start_pos = rgd.X.ncol - 2;
                    } else
                    {
                        start_pos = 0;
                    }
                    if (input_var.getInteraction() != 0
                            && !input_var.getAllcov()
                            && input_var.getNgpreds() != 2)
                    {
                        start_pos++;
                    }

                    if (model == 0)
                    {
                        end_pos = rgd.X.ncol;
                    } else
                    {
                        end_pos = rgd.X.ncol - 1;
                    }

                    if (input_var.getInteraction() != 0)
                    {
                        end_pos++;
                    }

                    for (int pos = start_pos; pos < end_pos; pos++)
                    {
                        *beta_sebeta[model] << input_var.getSep() << "nan"
                                << input_var.getSep() << "nan";
                    }

                    if (input_var.getNgpreds() == 2)
                    {
                        //Han Chen
#if !COXPH
                        if (!input_var.getAllcov()
                                && input_var.getInteraction() != 0)
                        {
                            if (model == 0)
                            {
                                *covvalue[model] << "nan" << input_var.getSep()
                                        << "nan";
                            } else
                            {
                                *covvalue[model] << "nan";
                            }
                        }
#endif
                        //Oct 26, 2009
                        *chi2[model] << "nan";
                    } else
                    { //ngpreds=1
                        if (input_var.getInverseFilename() == NULL)
                        {
                            //                     Han Chen
#if !COXPH
                            if (!input_var.getAllcov()
                                    && input_var.getInteraction() != 0)
                            {
                                *covvalue[model] << "nan";
                            }
#endif
                            //Oct 26, 2009
                            *chi2[model] << "nan";
                        }
                    }

                }
                //end of model cycle
                if (input_var.getNgpreds() == 2)
                {
                    //Han Chen
                    *outfile[0] << beta_sebeta[0]->str() << input_var.getSep();
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *outfile[0] << covvalue[0]->str() << input_var.getSep();
                    }
#endif
                    *outfile[0] << chi2[0]->str() << "\n";

                    *outfile[1] << beta_sebeta[1]->str() << input_var.getSep();
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *outfile[1] << covvalue[1]->str() << input_var.getSep();
                    }
#endif
                    *outfile[1] << chi2[1]->str() << "\n";

                    *outfile[2] << beta_sebeta[2]->str() << input_var.getSep();
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *outfile[2] << covvalue[2]->str() << input_var.getSep();
                    }
#endif
                    *outfile[2] << chi2[2]->str() << "\n";

                    *outfile[3] << beta_sebeta[3]->str() << input_var.getSep();
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *outfile[3] << covvalue[3]->str() << input_var.getSep();
                    }
#endif
                    *outfile[3] << chi2[3]->str() << "\n";

                    *outfile[4] << beta_sebeta[4]->str() << input_var.getSep();
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *outfile[4] << covvalue[4]->str() << input_var.getSep();
                    }
#endif
                    *outfile[4] << chi2[4]->str() << "\n";
                    //Oct 26, 2009

                } else  //ngpreds == 1
                {

                    if (input_var.getInverseFilename() == NULL)
                    {
                        //Han Chen
                        *outfile[0] << beta_sebeta[0]->str()
                                << input_var.getSep();
#if !COXPH
                        if (!input_var.getAllcov()
                                && input_var.getInteraction() != 0)
                        {
                            *outfile[0] << covvalue[0]->str()
                                    << input_var.getSep();
                        }
#endif
                        *outfile[0] << chi2[model]->str() << "\n";
                        //Oct 26, 2009
                    } else
                    {
                        *outfile[0] << beta_sebeta[0]->str() << "\n";
                    }
                }  // end ngpreds ==1
            }
        }
//
//End of All models output.
//
        else  //Only additive model. Only one output file
        {
            //Write mlinfo to output:
            int file = 0;
            write_mlinfo(outfile, file, mli, csnp, input_var, gcount, freq);
            int model = 0;
            if (poly)  //allel freq is not to rare
            {
#if LOGISTIC
                logistic_reg rd(rgd);
                if (input_var.getScore())
                rd.score(nrd.residuals, rgd, 0, CHOLTOL, model,
                        input_var.getInteraction(),
                        input_var.getNgpreds(),
                        invvarmatrix);
                else
                rd.estimate(rgd, 0, MAXITER, EPS, CHOLTOL, model,
                        input_var.getInteraction(),
                        input_var.getNgpreds(),
                        invvarmatrix,
                        input_var.getRobust());
#elif LINEAR
                //cout << (rgd.get_unmasked_data()).nids << " 1\n";
#if DEBUG
                rgd.X.print();
                rgd.Y.print();
#endif
                linear_reg rd(rgd);
#if DEBUG
                rgd.X.print();
                rgd.Y.print();
#endif
                //cout << (rgd.get_unmasked_data()).nids << " 2\n";
                if (input_var.getScore())
                {
#if DEBUG
                    cout << "input_var.getScore/n";
                    nrd.residuals.print();
                    cout << CHOLTOL << " <-CHOLTOL\n";
                    cout << model << " <-model\n";
                    cout << input_var.getInteraction()
                    << " <-input_var.getInteraction()\n";
                    cout << input_var.getNgpreds()
                    << " <-input_var.getNgpreds()\n";
                    invvarmatrix.print();
#endif
                    rd.score(nrd.residuals, rgd, 0, CHOLTOL, model,
                            input_var.getInteraction(), input_var.getNgpreds(),
                            invvarmatrix);
#if DEBUG
                    rd.beta.print();
                    cout << rd.chi2_score << " <-chi2_scoren\n";
                    rd.covariance.print();
                    rd.residuals.print();
                    rd.sebeta.print();
                    cout << rd.loglik << " <-logliken\n";
                    cout << rd.sigma2 << " <-sigma2\n";
#endif
                } else
                {
                    // if(input_var.getInverseFilename()== NULL)
                    // {
                    // cout << (rgd.get_unmasked_data()).nids << " 3\n";
#if DEBUG
                    cout << "rd.estimate\n";
                    cout << CHOLTOL << " <-CHOLTOL\n";
                    cout << model << " <-model\n";
                    cout << input_var.getInteraction()
                    << " <-input_var.getInteraction()\n";
                    cout << input_var.getNgpreds()
                    << " <-input_var.getNgpreds()\n";
                    cout << input_var.getRobust()
                    << " <-input_var.getRobust()\n";
                    cout << "start invarmatrix\n";
                    invvarmatrix.print();
                    cout << "end invarmatrix\n";
                    cout << rgd.is_interaction_excluded
                    << " <-rgd.is_interaction_excluded\n";
#endif
                    rd.estimate(rgd, 0, CHOLTOL, model,
                            input_var.getInteraction(), input_var.getNgpreds(),
                            invvarmatrix, input_var.getRobust());

#if DEBUG
                    cout << "rd.beta\n";
                    rd.beta.print();
                    cout << rd.chi2_score << " <-chi2_scoren\n";
                    cout << "rd.covariance\n";
                    rd.covariance.print();
                    cout << "rd.residuals\n";
                    rd.residuals.print();
                    cout << "rd.sebeta\n";
                    rd.sebeta.print();
                    cout << rd.loglik << " <-logliken\n";
                    cout << rd.sigma2 << " <-sigma2\n";
#endif
                    //cout << (rgd.get_unmasked_data()).nids << " 4\n";
                    //}
                    //else
                    //{
                    //   rd.mmscore(rgd, 0, CHOLTOL, model,
                    //              input_var.getInteraction(),
                    //              input_var.getNgpreds(), invvarmatrix);
                    //}
                }
#elif COXPH
                coxph_reg rd(rgd);
                rd.estimate(rgd, 0, MAXITER, EPS, CHOLTOL, model,
                        input_var.getInteraction(), true,
                        input_var.getNgpreds());
#endif

                if (!input_var.getAllcov() && input_var.getInteraction() == 0)
                {
                    start_pos = rd.beta.nrow - 1;
                } else if (!input_var.getAllcov()
                        && input_var.getInteraction() != 0)
                {
                    start_pos = rd.beta.nrow - 2;
                } else
                {
                    start_pos = 0;
                }
#if DEBUG
                cout << " start_pos;" << start_pos << "\n";
#endif
                for (int pos = start_pos; pos < rd.beta.nrow; pos++)
                {
                    *beta_sebeta[0] << input_var.getSep() << rd.beta[pos]
                            << input_var.getSep() << rd.sebeta[pos];
                    //Han Chen
#if !COXPH
                    if (input_var.getInverseFilename() == NULL
                            && !input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        if (pos > start_pos)
                        {
                            *covvalue[0] << rd.covariance[pos - 1];
                        }
                    }
#endif
                    //Oct 26, 2009
                }

                //calculate chi2
                //________________________________
                if (input_var.getInverseFilename() == NULL)
                {
#if DEBUG
                    cout << " inverse_filename == NULL" << "\n";
#endif
                    if (input_var.getScore() == 0)
                    {
                        *chi2[0] << rd.loglik; //2.*(rd.loglik-null_loglik);
                    } else
                    {
                        *chi2[0] << "nan"; //rd.chi2_score;
                    }
                }
                //________________________________
            } else //beta, sebeta = nan
            {
                if (!input_var.getAllcov() && input_var.getInteraction() == 0)
                    start_pos = rgd.X.ncol - 1;
                else if (!input_var.getAllcov()
                        && input_var.getInteraction() != 0)
                    start_pos = rgd.X.ncol - 2;
                else
                    start_pos = 0;

                end_pos = rgd.X.ncol;
                if (input_var.getInteraction() != 0)
                {
                    end_pos++;
                }
                if (input_var.getInteraction() != 0 && !input_var.getAllcov())
                {
                    start_pos++;
                }

                for (int pos = start_pos; pos < end_pos; pos++)
                {
                    *beta_sebeta[0] << input_var.getSep() << "nan"
                            << input_var.getSep() << "nan";
                }
                if (input_var.getInverseFilename() == NULL)
                {
                    //Han Chen
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        *covvalue[0] << "nan";
                    }
#endif
                    //Oct 26, 2009
                    *chi2[0] << "nan";
                }
            }

            if (input_var.getInverseFilename() == NULL)
            {
                //Han Chen
                *outfile[0] << beta_sebeta[0]->str() << input_var.getSep();
#if !COXPH
                if (!input_var.getAllcov() && input_var.getInteraction() != 0)
                {
                    *outfile[0] << covvalue[0]->str() << input_var.getSep();
                }
#endif
                *outfile[0] << chi2[model]->str() << "\n";
                //Oct 26, 2009
            } else
            {
                *outfile[0] << beta_sebeta[0]->str() << "\n";
#if DEBUG
                cout << "Se beta" << beta_sebeta[0] << "\n";
#endif
            }
        }
//
//End of Only additive model.
//
//clean chi2
        for (int i = 0; i < 5; i++)
        {
            beta_sebeta[i]->str("");
            //Han Chen
            covvalue[i]->str("");
            //Oct 26, 2009
            chi2[i]->str("");
        }
        update_progress_to_cmd_line(csnp, nsnps);
    }

    std::cout << "\b\b\b\b\b\b\b\b\b" << 100.;
    std::cout << "%... done\n";

//________________________________________________________________
//Maksim, 9 Jan, 2009

    for (unsigned int i = 0; i < outfile.size(); i++)
    {
        outfile[i]->close();
        delete outfile[i];
    }

//delete gtd;

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
