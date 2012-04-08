//=====================================================================================
//
//           Filename:  src/main.cpp
//
//        Description:  ProbABEL head file.
//
//            Version:  0.1-3
//            Created:  ---
//           Revision:  none
//
//             Author:  Yurii S. Aulchenko (cox, log, lin regressions)
//             Modified by: L.C. Karssen,
//                          Maksim V. Struchalin
// 
// modified 14-May-2009 by MVS:  interaction with SNP, interaction with SNP with exclusion of interacted covariates,
//                              mmscore implemented (poor me)
// modified 20-Jul-2009 by YSA: small changes, bug fix in mmscore option
// modified 22-Sep-2009 by YSA: "robust" option added
//
// Modified by Han Chen (hanchen@bu.edu) on Nov 9, 2009
// to extract the covariance between the estimate of beta(SNP) and the estimate of beta(interaction)
// based on src/main.cpp version 0.1-0 as of Oct 19, 2009
//
//            Company:  Department of Epidemiology, ErasmusMC Rotterdam, The Netherlands.
//              Email:  i.aoultchenko@erasmusmc.nl, m.struchalin@erasmusmc.nl
//
//=====================================================================================

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <sstream>

#include "mematrix.h"
#include "mematri1.h"
#include "data.h"
#include "reg1.h"
#include "comand_line_settings.h"

#define MAXITER 10
#define EPS 1.e-8
#define CHOLTOL 1.5e-12

int main(int argc, char * argv[])
{

    cmdvars input_var;
    input_var.set_variables(argc, argv);
    input_var.printinfo();

#if COXPH
    if (input_var.getScore())
    {
        fprintf(stderr,"\n\nOption --score is implemented for linear and logistic models only\n");
        exit(1);
    }
#endif
    //	if (allcov && ngpreds>1)
    //	{
    //		fprintf(stdout,"\n\nWARNING: --allcov allowed only for 1 predictor (MLDOSE)\n");
    //		allcov = 0;
    //	}

    mlinfo mli(input_var.getMlinfofilename(), input_var.getMapfilename());
    int nsnps = mli.nsnps;
    phedata phd;
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
        std::cerr << "error: Interaction parameter is out of range (ineraction="
                << input_var.getInteraction() << ") \n";
        exit(1);
    }

    //interaction--;

    //	if(input_var.getInverseFilename()!= NULL && phd.ncov > 1)
    //		{
    //		std::cerr<<"Error: In mmscore you can not use any covariates. You phenotype file must conatin id column and trait (residuals) only\n";
    //		exit(1);
    //		}

    //	if(input_var.getInverseFilename()!= NULL && (allcov == 1 || score == 1 || input_var.getInteraction()!= 0 || ngpreds==2))
    //		{
    //		std::cerr<<"Error: In mmscore you can use additive model without any inetractions only\n";
    //		exit(1);
    //		}

    mematrix<double> invvarmatrix;

    /*
     * now should be possible... delete this part later when everything works
     #if LOGISTIC
     if(input_var.getInverseFilename()!= NULL) {std::cerr<<"ERROR: mmscore is forbidden for logistic regression\n";exit(1);}
     #endif
     */

#if COXPH
    if(input_var.getInverseFilename()!= NULL)
    {
        std::cerr<<"ERROR: mmscore is forbidden for cox regression\n";
        exit(1);
    }
    if (input_var.getRobust())
    {
        std::cerr<<"ERROR: robust standard errors not implemented for Cox regression\n";
        exit(1);
    }
#endif

    if (input_var.getInverseFilename() != NULL)
    {
        std::cout << "you are running mmscore...\n";
    }

    std::cout << "Reading data ...";

    if (input_var.getInverseFilename() != NULL)
    {
        InvSigma inv(input_var.getInverseFilename(), &phd);
        invvarmatrix = inv.get_matrix();
        double par = 1.; //var(phd.Y)*phd.nids/(phd.nids-phd.ncov-1);
        invvarmatrix = invvarmatrix * par;
        std::cout << " loaded InvSigma ...";
        //	matrix.print();
    }

    std::cout.flush();

    gendata gtd;
    if (!input_var.getIsFvf())
        gtd.re_gendata(input_var.getGenfilename(), nsnps,
                input_var.getNgpreds(), phd.nids_all, phd.nids, phd.allmeasured,
                input_var.getSkipd(), phd.idnames);
    else
        gtd.re_gendata(input_var.getStrGenfilename(), nsnps,
                input_var.getNgpreds(), phd.nids_all, phd.nids, phd.allmeasured,
                phd.idnames);

    std::cout << " loaded genotypic data ...";

    /**
     if (input_var.getIsFvf())
     gendata gtd (str_genfilename,nsnps,input_var.getNgpreds(),phd.nids_all,phd.allmeasured,phd.idnames);
     else
     gendata gtd (input_var.getGenfilename(),nsnps,input_var.getNgpreds(),phd.nids_all,phd.nids,phd.allmeasured,skipd,phd.idnames);
     **/
    // estimate null model
    double null_loglik = 0.;
#if COXPH
    coxph_data nrgd=coxph_data(phd,gtd,-1,input_var.isIsInteractionExcluded());
#else
    regdata nrgd = regdata(phd, gtd, -1,input_var.isIsInteractionExcluded());
#endif

    std::cout << " loaded null data ...";

#if LOGISTIC
    logistic_reg nrd=logistic_reg(nrgd);
    nrd.estimate(nrgd,0,MAXITER,EPS,CHOLTOL,0,input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix, input_var.getRobust(), 1);
#elif LINEAR

    linear_reg nrd = linear_reg(nrgd);

    nrd.estimate(nrgd, 0, CHOLTOL, 0, input_var.getInteraction(),
            input_var.getNgpreds(), invvarmatrix, input_var.getRobust(), 1);
#elif COXPH
    coxph_reg nrd(nrgd);

    nrd.estimate(nrgd,0,MAXITER,EPS,CHOLTOL,0, input_var.getInteraction(), input_var.getNgpreds(), 1);
#endif
    null_loglik = nrd.loglik;

    std::cout << " estimated null model ...";

    // end null
#if COXPH
    coxph_data rgd(phd,gtd,0,input_var.isIsInteractionExcluded());
#else
    regdata rgd(phd, gtd, 0,input_var.isIsInteractionExcluded());
#endif

    std::cout << " formed regression object ...";

    std::cout << " done\n";
    std::cout.flush();

    //________________________________________________________________
    //Maksim, 9 Jan, 2009

    std::string outfilename_str(input_var.getOutfilename());
    std::vector<std::ofstream*> outfile;

    if (input_var.getNohead() != 1)
    {

        if (input_var.getNgpreds() == 2) //All models output. One file per each model
        {
            // open a file for output
            //_____________________

            for (int i = 0; i < 5; i++)
            {
                outfile.push_back(new std::ofstream());
            }

            outfile[0]->open((outfilename_str + "_2df.out.txt").c_str());
            outfile[1]->open((outfilename_str + "_add.out.txt").c_str());
            outfile[2]->open((outfilename_str + "_domin.out.txt").c_str());
            outfile[3]->open((outfilename_str + "_recess.out.txt").c_str());
            outfile[4]->open((outfilename_str + "_over_domin.out.txt").c_str());

            if (!outfile[0]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_2df.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[1]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_add.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[2]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_domin.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[3]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_recess.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[4]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_over_domin.out.txt" << "\n";
                exit(1);
            }
            //_____________________

            //Header
            //_____________________
            for (int i = 0; i < outfile.size(); i++)
            {
                (*outfile[i]) << "name" << input_var.getSep() << "A1"
                        << input_var.getSep() << "A2" << input_var.getSep()
                        << "Freq1" << input_var.getSep() << "MAF"
                        << input_var.getSep() << "Quality" << input_var.getSep()
                        << "Rsq" << input_var.getSep() << "n"
                        << input_var.getSep() << "Mean_predictor_allele";
                if (input_var.getChrom() != "-1")
                    (*outfile[i]) << input_var.getSep()
                            << "input_var.getChrom()";
                if (input_var.getMapfilename() != NULL)
                    (*outfile[i]) << input_var.getSep() << "position";
            }
            //_____________________

            if (input_var.getAllcov()) //All covariates in output
            {
                for (int file = 0; file < outfile.size(); file++)
                    for (int i = 0; i < phd.n_model_terms - 1; i++)
                        *outfile[file] << input_var.getSep() << "beta_"
                                << phd.model_terms[i] << input_var.getSep()
                                << "sebeta_" << phd.model_terms[i];
            }
            *outfile[0] << input_var.getSep() << "beta_SNP_A1A2"
                    << input_var.getSep() << "beta_SNP_A1A1"
                    << input_var.getSep() << "sebeta_SNP_A1A2"
                    << input_var.getSep() << "sebeta_SNP_A1A1";

            *outfile[1] << input_var.getSep() << "beta_SNP_addA1"
                    << input_var.getSep() << "sebeta_SNP_addA1";
            *outfile[2] << input_var.getSep() << "beta_SNP_domA1"
                    << input_var.getSep() << "sebeta_SNP_domA1";
            *outfile[3] << input_var.getSep() << "beta_SNP_recA1"
                    << input_var.getSep() << "sebeta_SNP_recA1";
            *outfile[4] << input_var.getSep() << "beta_SNP_odom"
                    << input_var.getSep() << "sebeta_SNP_odom";

            if (input_var.getInteraction() != 0)
            {
                //Han Chen
                *outfile[0] << input_var.getSep() << "beta_SNP_A1A2_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep() << "sebeta_SNP_A1A2_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep() << "beta_SNP_A1A1_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep() << "sebeta_SNP_A1A1_"
                        << phd.model_terms[interaction_cox];
#if !COXPH
                if (input_var.getInverseFilename() == NULL
                        && !input_var.getAllcov())
                {
                    *outfile[0] << input_var.getSep() << "cov_SNP_A1A2_int_SNP_"
                            << phd.model_terms[interaction_cox]
                            << input_var.getSep() << "cov_SNP_A1A1_int_SNP_"
                            << phd.model_terms[interaction_cox];
                }
#endif
                //Oct 26, 2009
                for (int file = 1; file < outfile.size(); file++)
                {
                    *outfile[file] << input_var.getSep() << "beta_SNP_"
                            << phd.model_terms[interaction_cox]
                            << input_var.getSep() << "sebeta_SNP_"
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
            *outfile[0] << input_var.getSep() << "loglik\n"; //"chi2_SNP_2df\n";
            *outfile[1] << input_var.getSep() << "loglik\n"; //"chi2_SNP_A1\n";
            *outfile[2] << input_var.getSep() << "loglik\n"; //"chi2_SNP_domA1\n";
            *outfile[3] << input_var.getSep() << "loglik\n"; //"chi2_SNP_recA1\n";
            *outfile[4] << input_var.getSep() << "loglik\n"; //"chi2_SNP_odom\n";

        }
        else //Only additive model. Only one output file
        {

            // open a file for output
            //_____________________
            //		if (outfilename != NULL)
            //	 		{
            outfile.push_back(
                    new std::ofstream(
                            (outfilename_str + "_add.out.txt").c_str()));
            //			}
            //		else
            //	 		{
            //			outfilename_str="regression_add.out.txt"; outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
            //			}

            if (!outfile[0]->is_open())
            {
                std::cerr << "Can not open file for writing: "
                        << outfilename_str << "\n";
                exit(1);
            }
            //_____________________

            //Header
            //_____________________
            *outfile[0] << "name" << input_var.getSep() << "A1"
                    << input_var.getSep() << "A2" << input_var.getSep()
                    << "Freq1" << input_var.getSep() << "MAF"
                    << input_var.getSep() << "Quality" << input_var.getSep()
                    << "Rsq" << input_var.getSep() << "n" << input_var.getSep()
                    << "Mean_predictor_allele";
            if (input_var.getChrom() != "-1")
                *outfile[0] << input_var.getSep() << "input_var.getChrom()";
            if (input_var.getMapfilename() != NULL)
                *outfile[0] << input_var.getSep() << "position";
            //_____________________

            if (input_var.getAllcov()) //All covariates in output
            {
                for (int i = 0; i < phd.n_model_terms - 1; i++)
                {
                    *outfile[0] << input_var.getSep() << "beta_"
                            << phd.model_terms[i] << input_var.getSep()
                            << "sebeta_" << phd.model_terms[i];
                }
                *outfile[0] << input_var.getSep() << "beta_SNP_add"
                        << input_var.getSep() << "sebeta_SNP_add";
            }
            else //Only beta, sebeta for additive model go to output file
            {
                *outfile[0] << input_var.getSep() << "beta_SNP_add"
                        << input_var.getSep() << "sebeta_SNP_add";
            }
            if (input_var.getInteraction() != 0)
            {
                *outfile[0] << input_var.getSep() << "beta_SNP_"
                        << phd.model_terms[interaction_cox]
                        << input_var.getSep() << "sebeta_SNP_"
                        << phd.model_terms[interaction_cox];
            }

            if (input_var.getInverseFilename() == NULL)
            //Han Chen
            {
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
    }
    else
    {
        if (input_var.getNgpreds() == 2) //All models output. One file per each model
        {
            // open a file for output
            //_____________________
            //		if (outfilename==NULL)
            //			{
            //			outfilename_str="regression";
            //			}

            for (int i = 0; i < 5; i++)
            {
                outfile.push_back(new std::ofstream());
            }

            outfile[0]->open((outfilename_str + "_2df.out.txt").c_str());
            outfile[1]->open((outfilename_str + "_add.out.txt").c_str());
            outfile[2]->open((outfilename_str + "_domin.out.txt").c_str());
            outfile[3]->open((outfilename_str + "_recess.out.txt").c_str());
            outfile[4]->open((outfilename_str + "_over_domin.out.txt").c_str());

            if (!outfile[0]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_2df.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[1]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_add.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[2]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_domin.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[3]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_recess.out.txt" << "\n";
                exit(1);
            }
            if (!outfile[4]->is_open())
            {
                std::cerr << "Cannot open file for writing: "
                        << outfilename_str + "_over_domin.out.txt" << "\n";
                exit(1);
            }
        }
        else
        {
            // open a file for output
            //_____________________
            //		if (outfilename != NULL)
            //	 		{
            outfile.push_back(
                    new std::ofstream(
                            (outfilename_str + "_add.out.txt").c_str()));
            //			}
            //		else
            //	 		{
            //			outfilename_str="regression_add.out.txt"; outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
            //			}

            if (!outfile[0]->is_open())
            {
                std::cerr << "Can not open file for writing: "
                        << outfilename_str << "\n";
                exit(1);
            }

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
     outfile << input_var.getSep() << "beta_SNP_A1A2" << input_var.getSep() << "beta_SNP_A1A1" << input_var.getSep()
     << "sebeta_SNP_A1A2" << input_var.getSep() << "sebeta_SNP_a1A1" << input_var.getSep() << "chi2_SNP_2df"
     << input_var.getSep() << "beta_SNP_addA1" << input_var.getSep() << "sebeta_SNP_addA1" << input_var.getSep() << "chi2_SNP_addA1"
     << input_var.getSep() << "beta_SNP_domA1" << input_var.getSep() << "sebeta_SNP_domA1" << input_var.getSep() << "chi2_SNP_domA1"
     << input_var.getSep() << "beta_SNP_recA1" << input_var.getSep() << "sebeta_SNP_recA1" << input_var.getSep() << "chi2_SNP_recA1"
     << input_var.getSep() << "beta_SNP_odom" << input_var.getSep() << "sebeta_SNP_odom" << input_var.getSep() << "chi2_SNP_odom\n";
     }
     else
     {
     outfile << input_var.getSep() << "beta_SNP_add" << input_var.getSep() << "sebeta_SNP_add" << input_var.getSep() << "chi2_SNP_add\n";
     }
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
            //		freq = ((gtd.G).column_mean(csnp*2)*2.+(gtd.G).column_mean(csnp*2+1))/2.;
            gtd.get_var(csnp * 2, snpdata1);
            gtd.get_var(csnp * 2 + 1, snpdata2);
            for (int ii = 0; ii < gtd.nids; ii++)
                if (!isnan(snpdata1[ii]) && !isnan(snpdata2[ii]))
                {
                    gcount++;
                    freq += snpdata1[ii] + snpdata2[ii] * 0.5;
                }
        }
        else
        {
            //		freq = (gtd.G).column_mean(csnp)/2.;
            gtd.get_var(csnp, snpdata1);
            for (int ii = 0; ii < gtd.nids; ii++)
                if (!isnan(snpdata1[ii]))
                {
                    gcount++;
                    freq += snpdata1[ii] * 0.5;
                }
        }
        freq /= (double) gcount;
        int poly = 1;
        if (fabs(freq) < 1.e-16 || fabs(1. - freq) < 1.e-16)
            poly = 0;
        if (fabs(mli.Rsq[csnp]) < 1.e-16)
            poly = 0;

        if (input_var.getNgpreds() == 2) //All models output. One file per each model
        {
            //Write mlinfo to output:
            for (int file = 0; file < outfile.size(); file++)
            {
                *outfile[file] << mli.name[csnp] << input_var.getSep()
                        << mli.A1[csnp] << input_var.getSep() << mli.A2[csnp]
                        << input_var.getSep() << mli.Freq1[csnp]
                        << input_var.getSep() << mli.MAF[csnp]
                        << input_var.getSep() << mli.Quality[csnp]
                        << input_var.getSep() << mli.Rsq[csnp]
                        << input_var.getSep() << gcount << input_var.getSep()
                        << freq;
                if (input_var.getChrom() != "-1")
                    *outfile[file] << input_var.getSep()
                            << input_var.getChrom();
                if (input_var.getMapfilename() != NULL)
                    *outfile[file] << input_var.getSep() << mli.map[csnp];
            }

            for (int model = 0; model < maxmod; model++)
            {
                if (poly) //allel freq is not to rare
                {
#if LOGISTIC
                    logistic_reg rd(rgd);
                    if (input_var.getScore())
                    rd.score(nrd.residuals,rgd,0,CHOLTOL,model,input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix);
                    else
                    rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model,input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix, input_var.getRobust());
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
                    rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model,input_var.getInteraction(), true, input_var.getNgpreds());
#endif

                    if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() == 0)
                        start_pos = rd.beta.nrow - 2;
                    else if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() != 0)
                        start_pos = rd.beta.nrow - 4;
                    else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() == 0)
                        start_pos = rd.beta.nrow - 1;
                    else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() != 0)
                        start_pos = rd.beta.nrow - 2;
                    else
                        start_pos = 0;

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
                                    if (pos > start_pos + 2)
                                    {
                                        *covvalue[model]
                                                << rd.covariance[pos - 3]
                                                << input_var.getSep()
                                                << rd.covariance[pos - 2];
                                    }
                                }
                                else
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
                    if (input_var.getScore() == 0)
                    {
                        //*chi2[model] << 2.*(rd.loglik-null_loglik);
                        *chi2[model] << rd.loglik;
                    }
                    else
                    {
                        //*chi2[model] << rd.chi2_score;
                        *chi2[model] << "nan";
                    }
                    //________________________________

                }
                else //beta, sebeta = nan
                {
                    if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() == 0)
                        start_pos = rgd.X.ncol - 2;
                    else if (!input_var.getAllcov() && model == 0
                            && input_var.getInteraction() != 0)
                        start_pos = rgd.X.ncol - 4;
                    else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() == 0)
                        start_pos = rgd.X.ncol - 1;
                    else if (!input_var.getAllcov() && model != 0
                            && input_var.getInteraction() != 0)
                        start_pos = rgd.X.ncol - 2;
                    else
                        start_pos = 0;

                    if (model == 0)
                    {
                        end_pos = rgd.X.ncol;
                    }
                    else
                    {
                        end_pos = rgd.X.ncol - 1;
                    }

                    if (input_var.getInteraction() != 0)
                        end_pos++;

                    for (int pos = start_pos; pos < end_pos; pos++)
                    {
                        *beta_sebeta[model] << input_var.getSep() << "nan"
                                << input_var.getSep() << "nan";
                    }
                    //Han Chen
#if !COXPH
                    if (!input_var.getAllcov()
                            && input_var.getInteraction() != 0)
                    {
                        if (model == 0)
                        {
                            *covvalue[model] << "nan" << input_var.getSep()
                                    << "nan";
                        }
                        else
                        {
                            *covvalue[model] << "nan";
                        }
                    }
#endif
                    //Oct 26, 2009
                    *chi2[model] << "nan";
                }
            } //end of moel cycle

            //Han Chen
            *outfile[0] << beta_sebeta[0]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[0] << covvalue[0]->str() << input_var.getSep();
            }
#endif
            *outfile[0] << chi2[0]->str() << "\n";
            *outfile[1] << beta_sebeta[1]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[1] << covvalue[1]->str() << input_var.getSep();
            }
#endif
            *outfile[1] << chi2[1]->str() << "\n";
            *outfile[2] << beta_sebeta[2]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[2] << covvalue[2]->str() << input_var.getSep();
            }
#endif
            *outfile[2] << chi2[2]->str() << "\n";
            *outfile[3] << beta_sebeta[3]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[3] << covvalue[3]->str() << input_var.getSep();
            }
#endif
            *outfile[3] << chi2[3]->str() << "\n";
            *outfile[4] << beta_sebeta[4]->str() << input_var.getSep();
#if !COXPH
            if (!input_var.getAllcov() && input_var.getInteraction() != 0)
            {
                *outfile[4] << covvalue[4]->str() << input_var.getSep();
            }
#endif
            *outfile[4] << chi2[4]->str() << "\n";
            //Oct 26, 2009

        }
        else //Only additive model. Only one output file
        {
            //Write mlinfo to output:
            *outfile[0] << mli.name[csnp] << input_var.getSep() << mli.A1[csnp]
                    << input_var.getSep() << mli.A2[csnp] << input_var.getSep();
            *outfile[0] << mli.Freq1[csnp] << input_var.getSep()
                    << mli.MAF[csnp] << input_var.getSep() << mli.Quality[csnp]
                    << input_var.getSep() << mli.Rsq[csnp]
                    << input_var.getSep();
            *outfile[0] << gcount << input_var.getSep() << freq;
            if (input_var.getChrom() != "-1")
                *outfile[0] << input_var.getSep() << input_var.getChrom();
            if (input_var.getMapfilename() != NULL)
                *outfile[0] << input_var.getSep() << mli.map[csnp];
            int model = 0;
            if (poly) //allel freq is not to rare
            {
#if LOGISTIC
                logistic_reg rd(rgd);
                if (input_var.getScore())
                rd.score(nrd.residuals,rgd,0,CHOLTOL,model, input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix);
                else
                rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix, input_var.getRobust());
#elif LINEAR
                //cout << (rgd.get_unmasked_data()).nids << " 1\n";
                linear_reg rd(rgd);
                //cout << (rgd.get_unmasked_data()).nids << " 2\n";
                if (input_var.getScore())
                    rd.score(nrd.residuals, rgd, 0, CHOLTOL, model,
                            input_var.getInteraction(), input_var.getNgpreds(),
                            invvarmatrix);
                else
                {
                    //					if(input_var.getInverseFilename()== NULL)
                    //						{
                    //cout << (rgd.get_unmasked_data()).nids << " 3\n";
                    rd.estimate(rgd, 0, CHOLTOL, model,
                            input_var.getInteraction(), input_var.getNgpreds(),
                            invvarmatrix, input_var.getRobust());
                    //cout << (rgd.get_unmasked_data()).nids << " 4\n";
                    //						}
                    //					else
                    //						{
                    //						rd.mmscore(rgd,0,CHOLTOL,model,  input_var.getInteraction(), input_var.getNgpreds(), invvarmatrix);
                    //						}
                }
#elif COXPH
                coxph_reg rd(rgd);
                rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, input_var.getInteraction(), true, input_var.getNgpreds());
#endif

                if (!input_var.getAllcov() && input_var.getInteraction() == 0)
                    start_pos = rd.beta.nrow - 1;
                else if (!input_var.getAllcov()
                        && input_var.getInteraction() != 0)
                    start_pos = rd.beta.nrow - 2;
                else
                    start_pos = 0;

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
                    if (input_var.getScore() == 0)
                    {
                        *chi2[0] << rd.loglik; //2.*(rd.loglik-null_loglik);
                    }
                    else
                    {
                        *chi2[0] << "nan"; //rd.chi2_score;
                    }
                }
                //________________________________
            }
            else //beta, sebeta = nan
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
            }
            else
            {
                *outfile[0] << beta_sebeta[0]->str() << "\n";
            }
        }

        //clean chi2
        for (int i = 0; i < 5; i++)
        {
            beta_sebeta[i]->str("");
            //Han Chen
            covvalue[i]->str("");
            //Oct 26, 2009
            chi2[i]->str("");
        }

        if (csnp % 1000 == 0)
        {
            if (csnp == 0)
            {
                fprintf(stdout, "Analysis: %6.2f ...",
                        100. * double(csnp) / double(nsnps));
            }
            else
            {
                fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b%6.2f ...",
                        100. * double(csnp) / double(nsnps));
            }
            std::cout.flush();
        }

    }

    fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b%6.2f", 100.);

    fprintf(stdout, " ... done\n");

    //________________________________________________________________
    //Maksim, 9 Jan, 2009

    for (int i = 0; i < outfile.size(); i++)
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
        it++;
    }
    it = covvalue.begin();
    while (it != covvalue.end())
    {
        delete *it;
        it++;
    }
    it = chi2.begin();
    while (it != chi2.end())
    {
        delete *it;
        it++;
    }

    return (0);
}
