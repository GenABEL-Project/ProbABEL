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
//             Modified by: Maksim V. Struchalin, 
// 
// modified 14-May-2009 by MVS: interaction with SNP, interaction with SNP with exclusion of interacted covariates,
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
#include <getopt.h>
#include <vector>
#include <sstream>


#include "version.h"
#include "mematrix.h"
#include "mematri1.h"
#include "data.h"
#include "reg1.h"
#include "usage.h"

#define MAXITER 10
#define EPS 1.e-8
#define CHOLTOL 1.5e-12


bool is_interaction_excluded = false; //Oh Holy Matrix, forgive me for this!

int main(int argc, char * argv [])
{

	int next_option;
	const char * const short_options = "p:i:d:m:n:c:o:s:t:g:a:erlh:b:vu";
	//b - interaction parameter
	// ADD --fv FLAG (FILEVECTOR), IN WHICH CASE USE ALTERNATIVE
	// CONSTRUCTOR FOR GENDATA
	const struct option long_options [] =
	{
			{"pheno",  1, NULL, 'p'},
			{"info",   1, NULL, 'i'},
			{"dose",   1, NULL, 'd'},
			{"map",	   1, NULL, 'm'},
			{"nids",   1, NULL, 'n'},
			{"chrom",  1, NULL, 'c'},
			{"out",    1, NULL, 'o'},
			{"skipd",  1, NULL, 's'},
			{"ntraits",1, NULL, 't'},
			{"ngpreds",1, NULL, 'g'},
			{"separat",1, NULL, 'a'},
			{"score",  0, NULL, 'r'},
			{"no-head",0, NULL, 'e'},
			{"allcov", 0, NULL, 'l'},
			{"help",   0, NULL, 'h'},
			{"interaction",   1, NULL, 'b'},
			{"interaction_only",   1, NULL, 'k'},
			{"mmscore",   1, NULL, 'v'},
			{"robust",   0, NULL, 'u'},
			{NULL  ,   0, NULL, 0  }
	};
	char * program_name = argv[0];

	char *phefilename = NULL;
	char *mlinfofilename = NULL;
	char *genfilename = NULL;
	char *mapfilename = NULL;
	char *outfilename = NULL;
	char *inverse_filename = NULL;
	string sep = " ";
	int nohead=0;
	int score=0;
	int npeople=-1;
	int ngpreds=1;
	int interaction=0;
	int interaction_excluded=0;
	int robust = 0;
	string chrom = "-1";
	int neco[] = {0,0,0};
	bool iscox=false;
	int isFVF = 0;
#if COXPH
	int noutcomes = 2;
	iscox=true;
#else
	int noutcomes = 1;
#endif
	int skipd = 2;
	int allcov = 0;
	do
	{
		next_option = getopt_long(argc,argv,short_options,long_options,NULL);
		switch (next_option)
		{
		case 'h': print_help(program_name,0);
		case 'p':
			phefilename = optarg;
			neco[0]=1;
			break;
		case 'i':
			mlinfofilename = optarg;
			neco[1]=1;
			break;
		case 'd':
			genfilename = optarg;
			neco[2]=1;
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
			sep = optarg;
			break;
		case 'e':
			nohead=1;
			break;
		case 'r':
			score=1;
			break;
		case 'l':
			allcov=1;
			break;
		case 'b':
			interaction=atoi(optarg);
			break;
		case 'k':
			interaction_excluded=atoi(optarg);
			break;
		case 'v':
			inverse_filename=optarg;
			break;
		case 'u':
			robust=1;
			break;


		case '?': print_usage(program_name,1);
		case -1 : break;
		default: abort();
		}
	}
	while (next_option != -1);


	fprintf(stdout,"ProbABEL v. %s (%s) (C) Yurii Aulchenko, Maksim Struchalin, EMCR\n\n",VERSION,DATE);
	if (neco[0]!=1 || neco[1]!=1 || neco[2]!=1)
	{
		print_usage(program_name,1);
	}

	if (score) {
		cout << "option --score supressed from v 0.1-6\n";
		exit(1);
	}

	string str_genfilename = genfilename;
	if (str_genfilename.find(FILEVECTOR_INDEX_FILE_SUFFIX)!=string::npos ||
			str_genfilename.find(FILEVECTOR_DATA_FILE_SUFFIX)!=string::npos) isFVF = 1;

	fprintf(stdout,"Options in effect:\n");
	fprintf(stdout,"\t --pheno   = %s\n",phefilename);
	fprintf(stdout,"\t --info    = %s\n",mlinfofilename);
	fprintf(stdout,"\t --dose    = %s\n",genfilename);
	if (isFVF) fprintf(stdout,"\t             (using FVF data)\n");
	fprintf(stdout,"\t --ntraits = %d\n",noutcomes);
	fprintf(stdout,"\t --ngpreds = %d\n",ngpreds);
	fprintf(stdout,"\t --interaction = %d\n",interaction);
	fprintf(stdout,"\t --interaction_only = %d\n",interaction_excluded);

	if (inverse_filename != NULL) fprintf(stdout,"\t --mmscore = %s\n",inverse_filename);
	else fprintf(stdout,"\t --mmscore     = not in output\n");
	//	fprintf(stdout,"\t --mmscore = %s\n",inverse_filename);

	if (mapfilename != NULL) fprintf(stdout,"\t --map     = %s\n",mapfilename);
	else fprintf(stdout,"\t --map     = not in output\n");
	if (npeople>0) fprintf(stdout,"\t --nids    = %d\n",npeople);
	else fprintf(stdout,"\t --nids    = estimated from data\n");
	if (chrom != "-1") cout << "\t --chrom   = " << chrom << "\n";
	else cout << "\t --chrom   = not in output\n";
	if (outfilename != NULL ) fprintf(stdout,"\t --out     = %s\n",outfilename);
	else fprintf(stdout,"\t --out     = regression.out.txt\n");
	fprintf(stdout,"\t --skipd   = %d\n",skipd);
	cout << "\t --separat = \"" << sep << "\"\n";
	if (score)
		fprintf(stdout,"\t --score   = ON\n");
	else
		fprintf(stdout,"\t --score   = OFF\n");
	if (nohead)
		fprintf(stdout,"\t --nohead  = ON\n");
	else
		fprintf(stdout,"\t --nohead  = OFF\n");
	if (allcov)
		fprintf(stdout,"\t --allcov  = ON\n");
	else
		fprintf(stdout,"\t --allcov  = OFF\n");
	if (robust)
		fprintf(stdout,"\t --robust  = ON\n");
	else
		fprintf(stdout,"\t --robust  = OFF\n");

	if (ngpreds!=1 && ngpreds!=2) 
	{
		fprintf(stderr,"\n\n--ngpreds should be 1 for MLDOSE or 2 for MLPROB\n");
		exit(1);
	}

	if(interaction_excluded != 0)
	{
		interaction = interaction_excluded; //ups
		is_interaction_excluded = true;
	}

#if COXPH
	if (score) 
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

	mlinfo mli(mlinfofilename,mapfilename);
	int nsnps = mli.nsnps;
	phedata phd(phefilename,noutcomes,npeople, interaction, iscox);

	int interaction_cox = interaction;
#if COXPH
	interaction_cox--;
#endif
	if(interaction < 0 || interaction > phd.ncov || interaction_cox > phd.ncov)
	{
		std::cerr << "error: Interaction parameter is out of range (ineraction="<<interaction<<") \n";
		exit(1);
	}

	//interaction--;

	//	if(inverse_filename != NULL && phd.ncov > 1)
	//		{
	//		std::cerr<<"Error: In mmscore you can not use any covariates. You phenotype file must conatin id column and trait (residuals) only\n";
	//		exit(1);
	//		}

	//	if(inverse_filename != NULL && (allcov == 1 || score == 1 || interaction != 0 || ngpreds==2))
	//		{
	//		std::cerr<<"Error: In mmscore you can use additive model without any inetractions only\n";
	//		exit(1);
	//		}




	mematrix<double> invvarmatrix;

#if LOGISTIC		
	if(inverse_filename != NULL) {std::cerr<<"ERROR: mmscore is forbidden for logistic regression\n";exit(1);}
#endif

#if COXPH
	if(inverse_filename != NULL) {std::cerr<<"ERROR: mmscore is forbidden for cox regression\n";exit(1);}
	if (robust) {std::cerr<<"ERROR: robust standard errors not implemented for Cox regression (drop us e-mail if you really need that)\n";exit(1);}
#endif


	if(inverse_filename != NULL)
	{
		std::cout<<"you are runing mmscore...\n";
	}

	std::cout << "Reading data ...";



	if(inverse_filename != NULL)
	{
		InvSigma inv(inverse_filename, &phd);
		invvarmatrix = inv.get_matrix();
		double par=1.; //var(phd.Y)*phd.nids/(phd.nids-phd.ncov-1);
		invvarmatrix = invvarmatrix*par;
		std::cout << " loaded InvSigma ...";
		//	matrix.print();
	}




	std::cout.flush();

	gendata gtd;
	if (!isFVF)
		gtd.re_gendata(genfilename,nsnps,ngpreds,phd.nids_all,phd.nids,phd.allmeasured,skipd,phd.idnames);
	else
		gtd.re_gendata(str_genfilename,nsnps,ngpreds,phd.nids_all,phd.nids,phd.allmeasured,phd.idnames);

	std::cout << " loaded genotypic data ...";

	/**
	if (isFVF)
		gendata gtd (str_genfilename,nsnps,ngpreds,phd.nids_all,phd.allmeasured,phd.idnames);
	else
		gendata gtd (genfilename,nsnps,ngpreds,phd.nids_all,phd.nids,phd.allmeasured,skipd,phd.idnames);
	 **/
	// estimate null model
	double null_loglik=0.;
#if COXPH
	coxph_data nrgd(phd,gtd,-1);
#else 
	regdata nrgd(phd,gtd,-1);
#endif

	std::cout << " loaded null data ...";

#if LOGISTIC
	logistic_reg nrd(nrgd);
	nrd.estimate(nrgd,0,MAXITER,EPS,CHOLTOL,0, interaction, ngpreds, invvarmatrix, robust, 1);
#elif LINEAR
	//fprintf(stdout,"HERE 0\n");
	linear_reg nrd(nrgd);
	//fprintf(stdout,"HERE 0\n");
	nrd.estimate(nrgd,0,CHOLTOL,0, interaction, ngpreds, invvarmatrix, robust, 1);
	//fprintf(stdout,"HERE 0\n");
#elif COXPH
	coxph_reg nrd(nrgd);
//	std::cout << " !!! AAA !!! ";
	nrd.estimate(nrgd,0,MAXITER,EPS,CHOLTOL,0, interaction, ngpreds, 1);
#endif
	null_loglik = nrd.loglik;

	std::cout << " estimated null model ...";

	// end null
#if COXPH
	coxph_data rgd(phd,gtd,0);
#else 
	regdata rgd(phd,gtd,0);
#endif

	std::cout << " formed regression object ...";

	std::cout << " done\n";
	std::cout.flush();


	//fprintf(stdout,"HERE 1\n");





	//________________________________________________________________
	//Maksim, 9 Jan, 2009

	if (outfilename==NULL)
	{
		outfilename="regression";
	}

	std::string outfilename_str(outfilename);
	std::vector<std::ofstream*> outfile;


	if (nohead!=1)
	{

		if(ngpreds==2) //All models output. One file per each model
		{
			// open a file for output
			//_____________________



			for(int i=0 ; i<5 ; i++)
			{
				outfile.push_back(new std::ofstream());
			}

			outfile[0]->open((outfilename_str+"_2df.out.txt").c_str());
			outfile[1]->open((outfilename_str+"_add.out.txt").c_str());
			outfile[2]->open((outfilename_str+"_domin.out.txt").c_str());
			outfile[3]->open((outfilename_str+"_recess.out.txt").c_str());
			outfile[4]->open((outfilename_str+"_over_domin.out.txt").c_str());



			if (!outfile[0]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_2df.out.txt" << "\n"; exit(1);}
			if (!outfile[1]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_add.out.txt" << "\n"; exit(1);}
			if (!outfile[2]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_domin.out.txt" << "\n"; exit(1);}
			if (!outfile[3]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_recess.out.txt" << "\n"; exit(1);}
			if (!outfile[4]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_over_domin.out.txt" << "\n"; exit(1);}
			//_____________________


			//Header
			//_____________________
			for(int i=0 ; i<outfile.size() ; i++)
			{
				(*outfile[i]) << "name" << sep << "A1" << sep << "A2" << sep << "Freq1" << sep << "MAF" << sep << "Quality" << sep << "Rsq"
						<< sep << "n" << sep << "Mean_predictor_allele";
				if (chrom != "-1") (*outfile[i]) << sep << "chrom";
				if (mapfilename != NULL) (*outfile[i]) << sep << "position";
			}
			//_____________________

			if(allcov) //All covariates in output
			{
				for (int file=0; file<outfile.size() ; file++)
					for (int i =0; i<phd.n_model_terms-1;i++)
						*outfile[file] << sep << "beta_" << phd.model_terms[i] << sep << "sebeta_" << phd.model_terms[i];
			}
			*outfile[0] << sep << "beta_SNP_A1A2" << sep << "beta_SNP_A1A1" << sep << "sebeta_SNP_A1A2" << sep << "sebeta_SNP_A1A1";
			*outfile[1] << sep << "beta_SNP_addA1" << sep << "sebeta_SNP_addA1";
			*outfile[2] << sep << "beta_SNP_domA1" << sep << "sebeta_SNP_domA1";
			*outfile[3] << sep << "beta_SNP_recA1" << sep << "sebeta_SNP_recA1";
			*outfile[4] << sep << "beta_SNP_odom" << sep << "sebeta_SNP_odom";

			if(interaction != 0)
			{
				//Han Chen
				*outfile[0] << sep << "beta_SNP_A1A2_" << phd.model_terms[interaction_cox] << sep << "sebeta_SNP_A1A2_" << phd.model_terms[interaction_cox]
				                                                                                                                           << sep << "beta_SNP_A1A1_" << phd.model_terms[interaction_cox] << sep << "sebeta_SNP_A1A1_" << phd.model_terms[interaction_cox];
#if !COXPH
				if(inverse_filename == NULL && !allcov) *outfile[0] << sep << "cov_SNP_A1A2_int_SNP_" << phd.model_terms[interaction_cox]
				                                                                                                         << sep << "cov_SNP_A1A1_int_SNP_" << phd.model_terms[interaction_cox];
#endif
				//Oct 26, 2009
				for (int file=1; file<outfile.size() ; file++)
				{		
					*outfile[file] << sep << "beta_SNP_" << phd.model_terms[interaction_cox]  << sep << "sebeta_SNP_" << phd.model_terms[interaction_cox];
					//Han Chen
#if !COXPH
					if(inverse_filename == NULL && !allcov) *outfile[file] << sep << "cov_SNP_int_SNP_" << phd.model_terms[interaction_cox];
#endif
					//Oct 26, 2009
				}
			}
			*outfile[0] << sep << "loglik\n"; //"chi2_SNP_2df\n";
			*outfile[1] << sep << "loglik\n"; //"chi2_SNP_A1\n";
			*outfile[2] << sep << "loglik\n"; //"chi2_SNP_domA1\n";
			*outfile[3] << sep << "loglik\n"; //"chi2_SNP_recA1\n";
			*outfile[4] << sep << "loglik\n"; //"chi2_SNP_odom\n";


		}
		else //Only additive model. Only one output file
		{

			// open a file for output
			//_____________________
			//		if (outfilename != NULL)
			//	 		{
			outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
			//			}
			//		else
			//	 		{
			//			outfilename_str="regression_add.out.txt"; outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
			//			}

			if (!outfile[0]->is_open())
			{
				std::cerr << "Can not open file for writing: " << outfilename_str << "\n";
				exit(1);
			}
			//_____________________

			//Header
			//_____________________
			*outfile[0] << "name" << sep << "A1" << sep << "A2" << sep << "Freq1" << sep << "MAF" << sep << "Quality" << sep << "Rsq" 
					<< sep << "n" << sep << "Mean_predictor_allele";
			if (chrom != "-1") *outfile[0] << sep << "chrom";
			if (mapfilename != NULL) *outfile[0] << sep << "position";
			//_____________________


			if(allcov) //All covariates in output
			{
				for (int i =0; i<phd.n_model_terms-1;i++)
					*outfile[0] << sep << "beta_" << phd.model_terms[i] << sep << "sebeta_" << phd.model_terms[i];

				*outfile[0] << sep << "beta_SNP_add" << sep << "sebeta_SNP_add";
			}
			else //Only beta, sebeta for additive model go to output file
			{
				*outfile[0] << sep << "beta_SNP_add" << sep << "sebeta_SNP_add";
			}
			if(interaction != 0) *outfile[0] << sep << "beta_SNP_" << phd.model_terms[interaction_cox] << sep << "sebeta_SNP_" << phd.model_terms[interaction_cox];

			if(inverse_filename == NULL)
				//Han Chen
			{
#if !COXPH
				if(interaction != 0 && !allcov) *outfile[0] << sep << "cov_SNP_int_SNP_" << phd.model_terms[interaction_cox];
#endif
				*outfile[0] << sep << "loglik"; //"chi2_SNP";
			}
			//Oct 26, 2009
			*outfile[0] << "\n";


		}
	}
	else
	{
		if(ngpreds==2) //All models output. One file per each model
		{
			// open a file for output
			//_____________________
			//		if (outfilename==NULL)
			//			{
			//			outfilename_str="regression";
			//			}



			for(int i=0 ; i<5 ; i++)
			{
				outfile.push_back(new std::ofstream());
			}

			outfile[0]->open((outfilename_str+"_2df.out.txt").c_str());
			outfile[1]->open((outfilename_str+"_add.out.txt").c_str());
			outfile[2]->open((outfilename_str+"_domin.out.txt").c_str());
			outfile[3]->open((outfilename_str+"_recess.out.txt").c_str());
			outfile[4]->open((outfilename_str+"_over_domin.out.txt").c_str());



			if (!outfile[0]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_2df.out.txt" << "\n"; exit(1);}
			if (!outfile[1]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_add.out.txt" << "\n"; exit(1);}
			if (!outfile[2]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_domin.out.txt" << "\n"; exit(1);}
			if (!outfile[3]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_recess.out.txt" << "\n"; exit(1);}
			if (!outfile[4]->is_open()) {std::cerr << "Can not open file for writing: " << outfilename_str+"_over_domin.out.txt" << "\n"; exit(1);}
		}
		else
		{
			// open a file for output
			//_____________________
			//		if (outfilename != NULL)
			//	 		{
			outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
			//			}
			//		else
			//	 		{
			//			outfilename_str="regression_add.out.txt"; outfile.push_back(new std::ofstream((outfilename_str+"_add.out.txt").c_str()));
			//			}

			if (!outfile[0]->is_open())
			{
				std::cerr << "Can not open file for writing: " << outfilename_str << "\n";
				exit(1);
			}

		}

	}

	//________________________________________________________________

	/*
	if (allcov)
		{
		if (score)
			{
			outfile << sep << "beta_mu"; // << sep << "beta_SNP_A1";
			outfile << sep << "sebeta_mu"; // << sep << "sebeta_SNP_A1";
			}
		else 
			{
			for (int i =0; i<phd.n_model_terms-1;i++) 
			outfile << sep << "beta_" << phd.model_terms[i] << sep << "sebeta_" << phd.model_terms[i];
			}
		if(interactio != 0) outfile << sep << "beta_SNP_" << phd.model_terms[interaction];
		}
	if (ngpreds==2) 
		{
		outfile << sep << "beta_SNP_A1A2" << sep << "beta_SNP_A1A1" << sep 
		<< "sebeta_SNP_A1A2" << sep << "sebeta_SNP_a1A1" << sep << "chi2_SNP_2df"
		<< sep << "beta_SNP_addA1" << sep << "sebeta_SNP_addA1" << sep << "chi2_SNP_addA1"
		<< sep << "beta_SNP_domA1" << sep << "sebeta_SNP_domA1" << sep << "chi2_SNP_domA1"
		<< sep << "beta_SNP_recA1" << sep << "sebeta_SNP_recA1" << sep << "chi2_SNP_recA1"
		<< sep << "beta_SNP_odom" << sep << "sebeta_SNP_odom" << sep << "chi2_SNP_odom\n";
		}	
	else 
		{
		outfile << sep << "beta_SNP_add" << sep << "sebeta_SNP_add" << sep << "chi2_SNP_add\n";
		}
	}
	 */
	//	exit(1);



	//________________________________________________________________
	//Maksim, 9 Jan, 2009



	int maxmod=5;
	int start_pos, end_pos;

	std::vector<std::ostringstream *> beta_sebeta;
	//Han Chen
	std::vector<std::ostringstream *> covvalue;
	//Oct 26, 2009
	std::vector<std::ostringstream *> chi2;

	for(int i=0 ; i<maxmod ; i++)
	{
		beta_sebeta.push_back(new std::ostringstream());
		//Han Chen
		covvalue.push_back(new std::ostringstream());
		//Oct 26, 2009
		chi2.push_back(new std::ostringstream());
	}



	for (int csnp=0;csnp<nsnps;csnp++)
	{

		rgd.update_snp(gtd,csnp);
		double freq = 0.;
		int gcount = 0;
		float snpdata1[gtd.nids];
		float snpdata2[gtd.nids];
		if (ngpreds==2)
		{
			//		freq = ((gtd.G).column_mean(csnp*2)*2.+(gtd.G).column_mean(csnp*2+1))/2.;
			gtd.get_var(csnp*2,snpdata1);
			gtd.get_var(csnp*2+1,snpdata2);
			for (int ii=0;ii<gtd.nids;ii++)
				if (!isnan(snpdata1[ii]) && !isnan(snpdata2[ii]))
				{gcount++;freq+=snpdata1[ii]+snpdata2[ii]*0.5;}
		} else {
			//		freq = (gtd.G).column_mean(csnp)/2.;
			gtd.get_var(csnp,snpdata1);
			for (int ii=0;ii<gtd.nids;ii++)
				if (!isnan(snpdata1[ii]))
				{gcount++;freq+=snpdata1[ii]*0.5;}
		}
		freq /= (double) gcount;
		int poly = 1;
		if (fabs(freq)<1.e-16 || fabs(1.-freq)<1.e-16) poly=0;
		if (fabs(mli.Rsq[csnp])<1.e-16) poly=0;



		if(ngpreds==2) //All models output. One file per each model
		{
			//Write mlinfo to output:
			for(int file=0 ; file<outfile.size() ; file++)
			{
				*outfile[file] << mli.name[csnp] << sep << mli.A1[csnp] << sep << mli.A2[csnp] << sep
						<< mli.Freq1[csnp] << sep << mli.MAF[csnp] << sep << mli.Quality[csnp] << sep
						<< mli.Rsq[csnp] << sep
						<< gcount << sep << freq;
				if (chrom != "-1") *outfile[file] << sep << chrom;
				if (mapfilename != NULL) *outfile[file] << sep << mli.map[csnp];
			}


			for(int model=0 ; model<maxmod ; model++)
			{		
				if(poly)//allel freq is not to rare
				{
#if LOGISTIC
					logistic_reg rd(rgd);
					if (score)
						rd.score(nrd.residuals,rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
					else
						rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, interaction, ngpreds, invvarmatrix, robust);
#elif LINEAR
					linear_reg rd(rgd);
					if(score)
						rd.score(nrd.residuals,rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
					else
					{
						//	rd.mmscore(rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
						rd.estimate(rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix, robust);
					}
#elif COXPH
					coxph_reg rd(rgd);
					rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, interaction, true, ngpreds);
#endif

					if(!allcov && model==0 && interaction==0) start_pos=rd.beta.nrow-2;
					else if(!allcov && model==0 && interaction!=0) start_pos=rd.beta.nrow-4;
					else if(!allcov && model!=0 && interaction==0) start_pos=rd.beta.nrow-1;
					else if(!allcov && model!=0 && interaction!=0) start_pos=rd.beta.nrow-2;
					else start_pos=0;


					for(int pos=start_pos ; pos<rd.beta.nrow ; pos++)
					{
						*beta_sebeta[model] << sep << rd.beta[pos] << sep << rd.sebeta[pos];
						//Han Chen
#if !COXPH
						if (inverse_filename == NULL && !allcov && interaction != 0)
						{
							if (pos>start_pos)
							{if (model==0)
							{if (pos>start_pos+2)
							{*covvalue[model] << rd.covariance[pos-3] << sep << rd.covariance[pos-2];}
							}
							else
							{*covvalue[model] << rd.covariance[pos-1];}
							}
						}
#endif
						//Oct 26, 2009
					}

					//calculate chi2
					//________________________________
					if (score==0)
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
					if(!allcov && model==0 && interaction==0) start_pos=rgd.X.ncol-2;
					else if(!allcov && model==0 && interaction!=0) start_pos=rgd.X.ncol-4;
					else if(!allcov && model!=0 && interaction==0) start_pos=rgd.X.ncol-1;
					else if(!allcov && model!=0 && interaction!=0) start_pos=rgd.X.ncol-2;
					else start_pos=0;

					if(model==0) {end_pos=rgd.X.ncol;}
					else {end_pos=rgd.X.ncol-1;}

					if(interaction!=0) end_pos++;

					for(int pos=start_pos ; pos<end_pos ; pos++)
					{
						*beta_sebeta[model] << sep << "nan" << sep << "nan";
					}
					//Han Chen
#if !COXPH
					if (!allcov && interaction !=0)
					{if (model==0)
					{*covvalue[model] << "nan" << sep << "nan";}
					else
					{*covvalue[model] << "nan";}
					}
#endif
					//Oct 26, 2009
					*chi2[model] << "nan";
				}
			}//end of moel cycle




			//Han Chen
			*outfile[0] << beta_sebeta[0]->str() << sep;
#if !COXPH
			if (!allcov && interaction !=0)
			{
				*outfile[0] << covvalue[0]->str() << sep;
			}
#endif
			*outfile[0] << chi2[0]->str() << "\n";
			*outfile[1] << beta_sebeta[1]->str() << sep;
#if !COXPH
			if (!allcov && interaction !=0)
			{
				*outfile[1] << covvalue[1]->str() << sep;
			}
#endif
			*outfile[1] << chi2[1]->str() << "\n";
			*outfile[2] << beta_sebeta[2]->str() << sep;
#if !COXPH
			if (!allcov && interaction !=0)
			{
				*outfile[2] << covvalue[2]->str() << sep;
			}
#endif
			*outfile[2] << chi2[2]->str() << "\n";
			*outfile[3] << beta_sebeta[3]->str() << sep;
#if !COXPH
			if (!allcov && interaction !=0)
			{
				*outfile[3] << covvalue[3]->str() << sep;
			}
#endif
			*outfile[3] << chi2[3]->str() << "\n";
			*outfile[4] << beta_sebeta[4]->str() << sep;
#if !COXPH
			if (!allcov && interaction !=0)
			{
				*outfile[4] << covvalue[4]->str() << sep;
			}
#endif
			*outfile[4] << chi2[4]->str() << "\n";
			//Oct 26, 2009





		}
		else //Only additive model. Only one output file
		{
			//Write mlinfo to output:
			*outfile[0] << mli.name[csnp] << sep << mli.A1[csnp] << sep << mli.A2[csnp] << sep;
			*outfile[0] << mli.Freq1[csnp] << sep << mli.MAF[csnp] << sep << mli.Quality[csnp] << sep << mli.Rsq[csnp] << sep;
			*outfile[0] << gcount << sep << freq;
			if (chrom != "-1") *outfile[0] << sep << chrom;
			if (mapfilename != NULL) *outfile[0] << sep << mli.map[csnp];
			int model=0;
			if(poly)//allel freq is not to rare
			{
#if LOGISTIC
				logistic_reg rd(rgd);
				if (score)
					rd.score(nrd.residuals,rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
				else
					rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, interaction, ngpreds, invvarmatrix, robust);
#elif LINEAR
				//cout << (rgd.get_unmasked_data()).nids << " 1\n";
				linear_reg rd(rgd);
				//cout << (rgd.get_unmasked_data()).nids << " 2\n";
				if (score)
					rd.score(nrd.residuals,rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
				else
				{
					//					if(inverse_filename == NULL)
					//						{
					//cout << (rgd.get_unmasked_data()).nids << " 3\n";
					rd.estimate(rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix, robust);
					//cout << (rgd.get_unmasked_data()).nids << " 4\n";
					//						}
					//					else
					//						{
					//						rd.mmscore(rgd,0,CHOLTOL,model, interaction, ngpreds, invvarmatrix);
					//						}
				}
#elif COXPH
				coxph_reg rd(rgd);
				rd.estimate(rgd,0,MAXITER,EPS,CHOLTOL,model, interaction, true, ngpreds);
#endif

				if(!allcov && interaction==0) start_pos=rd.beta.nrow-1;
				else if(!allcov && interaction!=0) start_pos=rd.beta.nrow-2;
				else start_pos=0;



				for(int pos=start_pos ; pos<rd.beta.nrow ; pos++)
				{
					*beta_sebeta[0] << sep << rd.beta[pos] << sep << rd.sebeta[pos];
					//Han Chen
#if !COXPH
					if (inverse_filename == NULL && !allcov && interaction != 0)
					{if (pos>start_pos)
					{*covvalue[0] << rd.covariance[pos-1];}
					}
#endif
					//Oct 26, 2009
				}





				//calculate chi2
				//________________________________
				if(inverse_filename == NULL)
				{
					if(score==0)
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
				if(!allcov && interaction==0) start_pos=rgd.X.ncol-1;
				else  if(!allcov && interaction!=0) start_pos=rgd.X.ncol-2;
				else start_pos=0;


				end_pos=rgd.X.ncol;
				if(interaction!=0) {end_pos++;}
				if(interaction!=0 && !allcov) {start_pos++;}



				for(int pos=start_pos ; pos<end_pos ; pos++)
				{
					*beta_sebeta[0] << sep << "nan" << sep << "nan";
				}
				if(inverse_filename == NULL)
				{
					//Han Chen
#if !COXPH
					if (!allcov && interaction !=0)
					{*covvalue[0] << "nan";}
#endif
					//Oct 26, 2009
					*chi2[0] << "nan";
				}
			}

			if(inverse_filename == NULL)
			{
				//Han Chen
				*outfile[0] << beta_sebeta[0]->str() << sep;
#if !COXPH
				if (!allcov && interaction !=0)
				{*outfile[0] << covvalue[0]->str() << sep;}
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
		for(int i=0 ; i<5 ; i++)
		{
			beta_sebeta[i]->str("");
			//Han Chen
			covvalue[i]->str("");
			//Oct 26, 2009
			chi2[i]->str("");
		}


		if (csnp % 1000 == 0)
		{
			if (csnp==0)
			{
				fprintf(stdout,"Analysis: %6.2f ...",100.*double(csnp)/double(nsnps));
			}
			else
			{
				fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b%6.2f ...",100.*double(csnp)/double(nsnps));
			}
			std::cout.flush();
		}

	}

	fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b%6.2f",100.);

	fprintf(stdout," ... done\n");

	//________________________________________________________________
	//Maksim, 9 Jan, 2009



	for(int i=0 ; i<outfile.size() ; i++)
	{
		outfile[i]->close();
		delete outfile[i];	
	}

	//delete gtd;

	return(0);
}
