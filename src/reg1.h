//==============================================================================
//
//           Filename:  src/reg1.h
//
//        Description:  ProbABEL
//
//            Version:  0.1-3
//            Created:  ---
//           Revision:  none
//  last modification:  11-Jan-2009
//
//             Author:  Yurii S. Aulchenko
//			  modified by: 	Maksim V. Struchalin, 11-Jan-2009
//
// Modified by Han Chen (hanchen@bu.edu) on Nov 9, 2009
// based on src/reg1.h version 0.2 as of Oct 19, 2009
//
//            Company:  ErasmusMC,
//                      Epidemiology & Biostatistics Department,
//                      Rotterdam,
//                      The Netherlands.
//              Email:  i.aoultchenko@erasmusmc.nl, m.struchalin@erasmusmc.nl
//
//==============================================================================

#include <cmath>
#include "cholesky.h"
#include "regdata.h"
#include "coxph_data.h"
extern "C" {
#include "survproto.h"
}

mematrix<double> apply_model(mematrix<double> &X, int model, int interaction,
		int ngpreds, bool is_interaction_excluded, bool iscox = false,
		int nullmodel = 0)
// model 0 = 2 df
// model 1 = additive 1 df
// model 2 = dominant 1 df
// model 3 = recessive 1 df
// model 4 = over-dominant 1 df

		{
#if DEBUG

	X.print();

	std::cout << "apply_model():model" << model << std::endl;
	std::cout << "apply_model():interaction" << interaction << std::endl;
	std::cout << "apply_model():ngpreds" << ngpreds << std::endl;
	std::cout << "apply_model():is_interaction_excluded"
	<< is_interaction_excluded << std::endl;
	std::cout << "apply_model():iscox" << iscox << std::endl;
	std::cout << "apply_model():nullmodel" << nullmodel << std::endl;
#endif
	if (model == 0) {
		if (interaction != 0 && !nullmodel) {
			if (ngpreds == 2) {
				mematrix<double> nX;
				nX.reinit(X.nrow, X.ncol + 2);
				int csnp_p1 = nX.ncol - 4;
				int csnp_p2 = nX.ncol - 3;
				int c1 = nX.ncol - 2;
				int c2 = nX.ncol - 1;
				for (int i = 0; i < X.nrow; i++)
					for (int j = 0; j < (X.ncol); j++)
						nX[i * nX.ncol + j] = X[i * X.ncol + j];

				for (int i = 0; i < nX.nrow; i++) {
					if (iscox) {
						//Maksim: interaction with SNP;;
						nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
								* X[i * X.ncol + interaction - 1];
						nX[i * nX.ncol + c2] = X[i * X.ncol + csnp_p2]
								* X[i * X.ncol + interaction - 1];
					} else {
						//Maksim: interaction with SNP;;
						nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
								* X[i * X.ncol + interaction];
						nX[i * nX.ncol + c2] = X[i * X.ncol + csnp_p2]
								* X[i * X.ncol + interaction];
					}
				}
				//________________________

				if (is_interaction_excluded) {
					mematrix<double> nX_without_interact_phe;
					nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);
					int col_new;
					for (int row = 0; row < nX.nrow; row++) {
						//Han Chen
						col_new = -1;
						for (int col = 0; col < nX.ncol; col++) {
							if (col != interaction && !iscox) {
								col_new++;
								nX_without_interact_phe[row
										* nX_without_interact_phe.ncol + col_new] =
										nX[row * nX.ncol + col];
							}
							if (col != interaction - 1 && iscox) {
								col_new++;
								nX_without_interact_phe[row
										* nX_without_interact_phe.ncol + col_new] =
										nX[row * nX.ncol + col];
							}
						} //interaction_only, model==0, ngpreds==2
						  //Oct 26, 2009
					}
					return nX_without_interact_phe;
				} //end of is_interaction_excluded
				  //________________________

				return (nX);
			}
			if (ngpreds == 1) {
				mematrix<double> nX;
				nX.reinit(X.nrow, X.ncol + 1);
				int csnp_p1 = nX.ncol - 2;
				int c1 = nX.ncol - 1;
				for (int i = 0; i < X.nrow; i++)
					for (int j = 0; j < (X.ncol); j++)
						nX[i * nX.ncol + j] = X[i * X.ncol + j];

				for (int i = 0; i < nX.nrow; i++) {
					if (iscox) {
						nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
								* X[i * X.ncol + interaction - 1]; //Maksim: interaction with SNP;;
					} else {
						nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
								* X[i * X.ncol + interaction]; //Maksim: interaction with SNP;;
					}
				}

				//________________________

				if (is_interaction_excluded) {
					int col_new;
					mematrix<double> nX_without_interact_phe;
					nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);
					for (int row = 0; row < nX.nrow; row++) {
						col_new = -1;
						for (int col = 0; col < nX.ncol; col++) {
							if (col != interaction && !iscox) {
								col_new++;
								nX_without_interact_phe[row
										* nX_without_interact_phe.ncol + col_new] =
										nX[row * nX.ncol + col];
							}
							if (col != interaction - 1 && iscox) {
								col_new++;
								nX_without_interact_phe[row
										* nX_without_interact_phe.ncol + col_new] =
										nX[row * nX.ncol + col];
							}

						}
					}
					return nX_without_interact_phe;
				} //end of is_interaction_excluded
				  //________________________

				return (nX);
			}
		} else {
			return (X);
		}
	}
	mematrix<double> nX;
	if (interaction != 0) {
		nX.reinit(X.nrow, (X.ncol));
	} else {
		nX.reinit(X.nrow, (X.ncol - 1));
	}
	int c1 = X.ncol - 2;
	int c2 = X.ncol - 1;
	for (int i = 0; i < X.nrow; i++)
		for (int j = 0; j < (X.ncol - 2); j++)
			nX[i * nX.ncol + j] = X[i * X.ncol + j];

	for (int i = 0; i < nX.nrow; i++) {
		if (model == 1)
			nX[i * nX.ncol + c1] = X[i * X.ncol + c1] + 2. * X[i * X.ncol + c2];
		else if (model == 2)
			nX[i * nX.ncol + c1] = X[i * X.ncol + c1] + X[i * X.ncol + c2];
		else if (model == 3)
			nX[i * nX.ncol + c1] = X[i * X.ncol + c2];
		else if (model == 4)
			nX[i * nX.ncol + c1] = X[i * X.ncol + c1];
		if (interaction != 0)
			nX[i * nX.ncol + c2] = X[i * nX.ncol + interaction]
					* nX[i * nX.ncol + c1]; //Maksim: interaction with SNP
	}
	//Han Chen

	if (is_interaction_excluded) {
		mematrix<double> nX_without_interact_phe;
		nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);
		int col_new;
		for (int row = 0; row < nX.nrow; row++) {
			col_new = -1;
			for (int col = 0; col < nX.ncol; col++) {
				if (col != interaction && !iscox) {
					col_new++;
					nX_without_interact_phe[row * nX_without_interact_phe.ncol
							+ col_new] = nX[row * nX.ncol + col];
				}
				if (col != interaction - 1 && iscox) {
					col_new++;
					nX_without_interact_phe[row * nX_without_interact_phe.ncol
							+ col_new] = nX[row * nX.ncol + col];
				}
			}
		}
		return nX_without_interact_phe;
	} //interaction_only, model!=0, ngpreds==2
	  //Oct 26, 2009
#if DEBUG
	nX.print();
#endif
	return nX;
}

mematrix<double> t_apply_model(mematrix<double> &X, int model, int interaction,
		int ngpreds, bool iscox, int nullmodel = 0) {
	mematrix<double> tmpX = transpose(X);
	mematrix<double> nX = apply_model(tmpX, model, interaction, ngpreds, iscox,
			nullmodel);
	mematrix<double> out = transpose(nX);
	return out;
}

class base_reg {
public:
	mematrix<double> beta;
	mematrix<double> sebeta;
	//Han Chen
	mematrix<double> covariance;
	//Oct 26, 2009
	mematrix<double> residuals;
	double sigma2;
	double loglik;
	double chi2_score;

	void base_score(mematrix<double> &resid, regdata &rdata, int verbose,
			double tol_chol, int model, int interaction, int ngpreds,
			masked_matrix invvarmatrix, int nullmodel = 0) {
		mematrix<double> oX = rdata.extract_genotypes();
		mematrix<double> X = apply_model(oX, model, interaction, ngpreds,
				rdata.is_interaction_excluded, false, nullmodel);
		beta.reinit(X.ncol, 1);
		sebeta.reinit(X.ncol, 1);
		double N = static_cast<double>(resid.nrow);

		mematrix<double> tX = transpose(X);
		if (invvarmatrix.length_of_mask != 0)
			tX = tX * invvarmatrix.masked_data;

		mematrix<double> u = tX * resid;
		mematrix<double> v = tX * X;
		mematrix<double> csum = column_sum(X);
		csum = transpose(csum) * csum;
		csum = csum * (1. / N);
		v = v - csum;
		// use cholesky to invert
		mematrix<double> v_i = v;
		cholesky2_mm(v_i, tol_chol);
		chinv2_mm(v_i);
		// before was
		// mematrix<double> v_i = invert(v);
		beta = v_i * u;
		double sr = 0.;
		double srr = 0.;
		for (int i = 0; i < resid.nrow; i++) {
			sr += resid[i];
			srr += resid[i] * resid[i];
		}
		double mean_r = sr / N;
		double sigma2_internal = (srr - N * mean_r * mean_r) / (N - beta.nrow);
		for (int i = 0; i < beta.nrow; i++)
			sebeta[i] = sqrt(v_i.get(i, i) * sigma2_internal);
		mematrix<double> chi2 = transpose(u) * v_i * u;
		chi2 = chi2 * (1. / sigma2_internal);
		chi2_score = chi2[0];
	}
};

class linear_reg: public base_reg {
public:

	linear_reg(regdata &rdatain) {
		regdata rdata = rdatain.get_unmasked_data();
		//fprintf(stdout,"linear_reg: %i %i %i\n",rdata.nids,(rdata.X).ncol,(rdata.Y).ncol);
		int length_beta = (rdata.X).ncol;
		beta.reinit(length_beta, 1);
		sebeta.reinit(length_beta, 1);
		//Han Chen
		if (length_beta > 1) {
			covariance.reinit(length_beta - 1, 1);
		}
		//Oct 26, 2009
		residuals.reinit(rdata.nids, 1);
		sigma2 = -1.;
		loglik = -9.999e+32;
		chi2_score = -1.;
	}
	~linear_reg() {
		//		delete beta;
		//		delete sebeta;
		//		delete residuals;
	}

//    mematrix<double> fill_invvar_matrix(mematrix<double>& invvarmatrix,
//            regdata& rdata, regdata& rdatain, mematrix<double> invvarmatrixin)
//    {
//        std::cout<<"invvarmatrixin"<<invvarmatrixin.ncol<<" "<<invvarmatrixin.ncol<<"\n";
//        std::cout<<"rdata.nids"<<rdata.nids<<"\n";
//
//        invvarmatrix.reinit(rdata.nids, rdata.nids);
//        int i1 = 0, j1 = 0;
//        for (int i = 0; i < rdata.nids; i++)
//            if (rdatain.masked_data[i] == 0)
//            {
//                j1 = 0;
//                for (int j = 0; j < rdata.nids; j++)
//                    if (rdatain.masked_data[j] == 0)
//                    {
//                        invvarmatrix.put(invvarmatrixin.get(i, j), i1, j1);
//                        j1++;
//                    }
//                i1++;
//            }
//
//        return invvarmatrix;
//    }

	void estimate(regdata& rdatain, int verbose, double tol_chol, int model,
			int interaction, int ngpreds, masked_matrix invvarmatrixin,
			int robust, int nullmodel = 0) {
		//suda ineraction parameter
		// model should come here
		regdata rdata = rdatain.get_unmasked_data();

		if (invvarmatrixin.length_of_mask != 0) {
			invvarmatrixin.update_mask(rdatain.masked_data);
			//  invvarmatrixin.masked_data->print();
		}
		if (verbose) {
			cout << rdata.is_interaction_excluded
					<< " <-irdata.is_interaction_excluded\n";
			printf("invvarmatrix:\n");
			invvarmatrixin.masked_data->print();
			printf("rdata.X:\n");
			rdata.X.print();
		}

		//fprintf(stdout,"estimate: %i %i %i %i\n",rdata.nids,(rdata.X).nrow,(rdata.X).ncol,(rdata.Y).ncol);
		mematrix<double> X = apply_model(rdata.X, model, interaction, ngpreds,
				rdata.is_interaction_excluded, false, nullmodel);

		if (verbose) {
			printf("X:\n");
			X.print();
		}
		int length_beta = X.ncol;
		beta.reinit(length_beta, 1);
		sebeta.reinit(length_beta, 1);
		//Han Chen
		if (length_beta > 1) {
			if (model == 0 && interaction != 0 && ngpreds == 2
					&& length_beta > 2) {
				covariance.reinit(length_beta - 2, 1);
			} else {
				covariance.reinit(length_beta - 1, 1);
			}
		}
		//Oct 26, 2009
		mematrix<double> tX = transpose(X);
		if (invvarmatrixin.length_of_mask != 0) {

			tX = tX * invvarmatrixin.masked_data;
			//!check if quicker
			//tX = productXbySymM(tX,invvarmatrix);
			// = invvarmatrix*X; std::cout<<"new tX.nrow="<<X.nrow<<" tX.ncol="<<X.ncol<<"\n";
		}

		mematrix<double> tXX = tX * X;

		//		double N = tXX.get(0,0);
		double N = X.nrow;

		//
		// use cholesky to invert
		//
		mematrix<double> tXX_i = tXX;
		cholesky2_mm(tXX_i, tol_chol);
		chinv2_mm(tXX_i);
		// before was
		// mematrix<double> tXX_i = invert(tXX);
		mematrix<double> tXY = tX * (rdata.Y);

		beta = tXX_i * tXY;
		if (verbose) {
			printf("tX:\n");
			tX.print();
			printf("tXX:\n");
			tXX.print();
			printf("chole tXX:\n");
			tXX_i.print();
			printf("tXX-1:\n");
			tXX_i.print();
			printf("tXY:\n");
			tXY.print();
			printf("beta:\n");
			(beta).print();
		}

		// now compute residual variance
		sigma2 = 0.;

		mematrix<double> ttX = transpose(tX);
		mematrix<double> sigma2_matrix = rdata.Y;
		mematrix<double> sigma2_matrix1 = ttX * beta;
//        printf("sigma2_matrix");
//        sigma2_matrix.print();
//
//        printf("sigma2_matrix1");
//        sigma2_matrix1.print();

		sigma2_matrix = sigma2_matrix - sigma2_matrix1;
//        printf("sigma2_matrix");
//        sigma2_matrix.print();
		static double val;

		//	std::cout<<"sigma2_matrix.nrow="<<sigma2_matrix.nrow<<"sigma2_matrix.ncol"<<sigma2_matrix.ncol<<"\n";

		for (int i = 0; i < sigma2_matrix.nrow; i++) {
			val = sigma2_matrix.get(i, 0);
//            printf("val = %f\n", val);
			sigma2 += val * val;
//            printf("sigma2+= = %f\n", sigma2);
		}

		double sigma2_internal = sigma2
				/ (N - static_cast<double>(length_beta));

		// now compute residual variance
		//		sigma2 = 0.;
		//		for (int i =0;i<(rdata.Y).nrow;i++)
		//			sigma2 += ((rdata.Y).get(i,0))*((rdata.Y).get(i,0));
		//		for (int i=0;i<length_beta;i++)
		//			sigma2 -= 2. * (beta.get(i,0)) * tXY.get(i,0);
		//		for (int i=0;i<(length_beta);i++)
		//		for (int j=0;j<(length_beta);j++)
		//			sigma2 += (beta.get(i,0)) * (beta.get(j,0)) * tXX.get(i,j);

		//	std::cout<<"sigma2="<<sigma2<<"\n";
		//	std::cout<<"sigma2_internal="<<sigma2_internal<<"\n";

		//		replaced for ML
		//		sigma2_internal	= sigma2/(N - double(length_beta) - 1);
//        printf("sigma2/=N = %f\n", sigma2);
		sigma2 /= N;

		//	std::cout<<"N="<<N<<", length_beta="<<length_beta<<"\n";

		if (verbose) {
			printf("sigma2 = %f\n", sigma2);
		}

		/*
		 loglik = 0.;
		 double ss=0;
		 for (int i=0;i<rdata.nids;i++) {
		 double resid = rdata.Y[i] - beta.get(0,0); // intercept
		 for (int j=1;j<beta.nrow;j++) resid -= beta.get(j,0)*X.get(i,j);
		 //	residuals[i] = resid;
		 ss += resid*resid;
		 }
		 sigma2 = ss/N;
		 */

		//cout << "estimate " << rdata.nids << "\n";
		//(rdata.X).print();
		//for (int i=0;i<rdata.nids;i++) cout << rdata.masked_data[i] << " ";
		//cout << endl;
		loglik = 0.;
		double halfrecsig2 = .5 / sigma2;
		for (int i = 0; i < rdata.nids; i++) {
			double resid = rdata.Y[i] - beta.get(0, 0); // intercept
			for (int j = 1; j < beta.nrow; j++)
				resid -= beta.get(j, 0) * X.get(i, j);
			residuals[i] = resid;
			loglik -= halfrecsig2 * resid * resid;
		}
		loglik -= static_cast<double>(rdata.nids) * log(sqrt(sigma2));

		//cout << "estimate " << rdata.nids << "\n";
		//
		//		ugly fix to the fact that if we do mmscore, sigma2 is already in the matrix...
		//		YSA, 2009.07.20
		//
		//cout << "estimate 0\n";

		if (invvarmatrixin.length_of_mask != 0)
			sigma2_internal = 1.0;

		mematrix<double> robust_sigma2(X.ncol, X.ncol);
		if (robust) {
			mematrix<double> XbyR = X;
			for (int i = 0; i < X.nrow; i++)
				for (int j = 0; j < X.ncol; j++) {
					double tmpval = XbyR.get(i, j) * residuals[i];
					XbyR.put(tmpval, i, j);
				}
			XbyR = transpose(XbyR) * XbyR;
			robust_sigma2 = tXX_i * XbyR;
			robust_sigma2 = robust_sigma2 * tXX_i;
		}

		//cout << "estimate 0\n";

		for (int i = 0; i < (length_beta); i++) {
			if (robust) {
//                cout << "estimate :robust\n";
				double value = sqrt(robust_sigma2.get(i, i));
				sebeta.put(value, i, 0);
				//Han Chen
				if (i > 0) {
					if (model == 0 && interaction != 0 && ngpreds == 2
							&& length_beta > 2) {
						if (i > 1) {
							double covval = robust_sigma2.get(i, i - 2);
							covariance.put(covval, i - 2, 0);
						}
					} else {
						double covval = robust_sigma2.get(i, i - 1);
						covariance.put(covval, i - 1, 0);
					}
				}
				//Oct 26, 2009
			} else {
//                cout << "estimate :non-robust\n";
				double value = sqrt(sigma2_internal * tXX_i.get(i, i));
				sebeta.put(value, i, 0);
				//Han Chen
				if (i > 0) {
					if (model == 0 && interaction != 0 && ngpreds == 2
							&& length_beta > 2) {
						if (i > 1) {
							double covval = sigma2_internal
									* tXX_i.get(i, i - 2);
							covariance.put(covval, i - 2, 0);
						}
					} else {
						double covval = sigma2_internal * tXX_i.get(i, i - 1);
						covariance.put(covval, i - 1, 0);
					}
				}
				//Oct 26, 2009
			}
		}
		//cout << "estimate E\n";
		if (verbose) {
			printf("sebeta (%d):\n", sebeta.nrow);
			sebeta.print();
		}
	}

	void score(mematrix<double> &resid, regdata &rdatain, int verbose,
			double tol_chol, int model, int interaction, int ngpreds,
			masked_matrix invvarmatrix, int nullmodel = 0) {
		regdata rdata = rdatain.get_unmasked_data();
		base_score(resid, rdata, verbose, tol_chol, model, interaction, ngpreds,
				invvarmatrix, nullmodel = 0);

	}
};

class logistic_reg: public base_reg {
public:
	int niter;

	logistic_reg(regdata &rdatain) {
		regdata rdata = rdatain.get_unmasked_data();
		int length_beta = (rdata.X).ncol;
		beta.reinit(length_beta, 1);
		sebeta.reinit(length_beta, 1);
		//Han Chen
		if (length_beta > 1) {
			covariance.reinit(length_beta - 1, 1);
		}
		//Oct 26, 2009
		residuals.reinit((rdata.X).nrow, 1);
		sigma2 = -1.;
		loglik = -9.999e+32; // should actually be MAX of the corresponding type
		niter = -1;
		chi2_score = -1.;
	}
	~logistic_reg() {
		//		delete beta;
		//		delete sebeta;
	}

	void estimate(regdata &rdatain, int verbose, int maxiter, double eps,
			double tol_chol, int model, int interaction, int ngpreds,
			masked_matrix invvarmatrixin, int robust, int nullmodel = 0) {
		// on the contrast to the 'linear' 'invvarmatrix' contains
		// inverse of correlation matrix (not the inverse of var-cov matrix)
		// h2.object$InvSigma * h.object2$h2an$estimate[length(h2$h2an$estimate)]
		// the inverse of var-cov matrix scaled by total variance
		regdata rdata = rdatain.get_unmasked_data();

		// a lot of code duplicated between linear and logistic...
		// e.g. a piece below...
		mematrix<double> invvarmatrix;
		if (invvarmatrixin.length_of_mask != 0) {
			invvarmatrixin.update_mask(rdatain.masked_data);

		}

		mematrix<double> X = apply_model(rdata.X, model, interaction, ngpreds,
				rdata.is_interaction_excluded, false, nullmodel);
		int length_beta = X.ncol;
		beta.reinit(length_beta, 1);
		sebeta.reinit(length_beta, 1);
		//Han Chen
		if (length_beta > 1) {
			if (model == 0 && interaction != 0 && ngpreds == 2
					&& length_beta > 2) {
				covariance.reinit(length_beta - 2, 1);
			} else {
				covariance.reinit(length_beta - 1, 1);
			}
		}
		//Oct 26, 2009
		mematrix<double> W((X).nrow, 1);
		mematrix<double> z((X).nrow, 1);
		mematrix<double> tXWX(length_beta, length_beta);
		mematrix<double> tXWX_i(length_beta, length_beta);
		mematrix<double> tXWz(length_beta, 1);

		double prev = (rdata.Y).column_mean(0);
		if (prev >= 1. || prev <= 0.) {
			fprintf(stderr, "prevalence not within (0,1)\n");
			exit(1);
		}
		for (int i = 0; i < length_beta; i++)
			beta.put(0., i, 0);
		beta.put(log(prev / (1. - prev)), 0, 0);

		mematrix<double> tX = transpose(X);

		if (invvarmatrix.nrow != 0 && invvarmatrix.ncol != 0) {
			tX = tX * invvarmatrix;
		}
		/*
		 fprintf(stdout,"\n");
		 fprintf(stdout,"X %f %f %f\n",X.get(0,0),X.get(0,1),X.get(0,2));
		 if (X.ncol==4) fprintf(stdout,"X[4] %f\n",X.get(0,3));
		 fprintf(stdout,"Inv %f %f %f\n",invvarmatrix.get(0,0),invvarmatrix.get(0,1),invvarmatrix.get(0,2));
		 if (X.ncol==4) fprintf(stdout,"X[4] %f\n",invvarmatrix.get(0,3));
		 fprintf(stdout,"tXInv %f %f %f\n",tX.get(0,0),tX.get(1,0),tX.get(2,0));
		 if (X.ncol==4) fprintf(stdout,"X[4] %f\n",tX.get(3,0));
		 */
		//TODO(maarten): remove this unused variable if there is not a reason to keep it
		//double N;
		niter = 0;
		double delta = 1.;
		double prevlik = 0.;
		while (niter < maxiter && delta > eps) {
			mematrix<double> eMu = (X) * beta;
			mematrix<double> eMu_us = eMu;
			for (int i = 0; i < eMu.nrow; i++) {
				double emu = eMu.get(i, 0);
				double value = emu;
				double zval = 0.;
				value = exp(value) / (1. + exp(value));
				residuals[i] = (rdata.Y).get(i, 0) - value;
				eMu.put(value, i, 0);
				W.put(value * (1. - value), i, 0);
				zval = emu
						+ (1. / (value * (1. - value)))
								* (((rdata.Y).get(i, 0)) - value);
				z.put(zval, i, 0);
			}

			mematrix<double> tmp = productMatrDiag(tX, W);
			if (verbose) {
				printf("tXW:\n");
				tmp.print();
			}
			mematrix<double> tXWX = tmp * (X);
			//N = tXWX.get(0, 0);

			if (verbose) {
				printf("tXWX:\n");
				tXWX.print();
			}
			//printf("tXWX:\n");tXWX.print();
			//
			// use cholesky to invert
			//
			tXWX_i = tXWX;
			//cholesky2_mm(tXWX_i,tol_chol);
			//if (verbose) {printf("chole tXWX:\n");tXWX_i.print();}
			//printf("chole tXWX:\n");tXWX_i.print();
			//chinv2_mm(tXWX_i);
			// was before
			tXWX_i = invert(tXWX);
			if (verbose) {
				printf("tXWX-1:\n");
				tXWX_i.print();
			}
			//fprintf(stdout,"*** tXWX_i\n");tXWX_i.print();
			mematrix<double> tmp1 = productMatrDiag(tX, W);
			mematrix<double> tXWz = tmp1 * z;
			if (verbose) {
				printf("tXWz:\n");
				tXWz.print();
			}
			beta = tXWX_i * tXWz;
			//fprintf(stdout,"*** res: %f %f %f\n",residuals[0],residuals[1],residuals[2]);
			//mematrix<double> txres = tx * residuals;
			//fprintf(stdout,"*** txres\n");txres.print();
			//beta = txwx_i* txres;
			if (verbose) {
				printf("beta:\n");
				beta.print();
			}
			//printf("beta:\n");beta.print();
			// compute likelihood
			prevlik = loglik;
			loglik = 0.;
			for (int i = 0; i < eMu.nrow; i++)
				loglik += rdata.Y[i] * eMu_us[i] - log(1. + exp(eMu_us[i]));
			delta = fabs(1. - (prevlik / loglik));
			niter++;
		}
		sigma2 = 0.;

		mematrix<double> robust_sigma2(X.ncol, X.ncol);
		if (robust) {
			mematrix<double> XbyR = X;
			for (int i = 0; i < X.nrow; i++)
				for (int j = 0; j < X.ncol; j++) {
					double tmpval = XbyR.get(i, j) * residuals[i];
					XbyR.put(tmpval, i, j);
				}
			XbyR = transpose(XbyR) * XbyR;
			robust_sigma2 = tXWX_i * XbyR;
			robust_sigma2 = robust_sigma2 * tXWX_i;
		}

		for (int i = 0; i < (length_beta); i++) {
			if (robust) {
				double value = sqrt(robust_sigma2.get(i, i));
				sebeta.put(value, i, 0);
				//Han Chen
				if (i > 0) {
					if (model == 0 && interaction != 0 && ngpreds == 2
							&& length_beta > 2) {
						if (i > 1) {
							double covval = robust_sigma2.get(i, i - 2);
							covariance.put(covval, i - 2, 0);
						}
					} else {
						double covval = robust_sigma2.get(i, i - 1);
						covariance.put(covval, i - 1, 0);
					}
				}
				//Oct 26, 2009
			} else {
				double value = sqrt(tXWX_i.get(i, i));
				sebeta.put(value, i, 0);
				//Han Chen
				if (i > 0) {
					if (model == 0 && interaction != 0 && ngpreds == 2
							&& length_beta > 2) {
						if (i > 1) {
							double covval = tXWX_i.get(i, i - 2);
							covariance.put(covval, i - 2, 0);
						}
					} else {
						double covval = tXWX_i.get(i, i - 1);
						covariance.put(covval, i - 1, 0);
					}
				}
				//Oct 26, 2009
			}
		}
		if (verbose) {
			printf("sebeta (%d):\n", sebeta.nrow);
			sebeta.print();
		}
		//printf("sebeta (%d):\n",beta.nrow);beta.print();
		//printf("sebeta (%d):\n",sebeta.nrow);sebeta.print();
		//exit(1);
	}
	// just a stupid copy from linear_reg
	void score(mematrix<double> &resid, regdata &rdata, int verbose,
			double tol_chol, int model, int interaction, int ngpreds,
			masked_matrix invvarmatrix, int nullmodel = 0) {
		base_score(resid, rdata, verbose, tol_chol, model, interaction, ngpreds,
				invvarmatrix, nullmodel = 0);
	}
};

//class coxph_reg
//{
//public:
//    mematrix<double> beta;
//    mematrix<double> sebeta;
//    mematrix<double> residuals;
//    double sigma2;
//    double loglik;
//    double chi2_score;
//    int niter;
//
//    coxph_reg(coxph_data &cdatain)
//    {
//        coxph_data cdata = cdatain.get_unmasked_data();
//        beta.reinit(cdata.X.nrow, 1);
//        sebeta.reinit(cdata.X.nrow, 1);
//        loglik = -9.999e+32;
//        sigma2 = -1.;
//        chi2_score = -1.;
//        niter = 0;
//    }
//    ~coxph_reg()
//    {
//        //		delete beta;
//        //		delete sebeta;
//    }
//    void estimate(coxph_data &cdatain, int verbose, int maxiter, double eps,
//            double tol_chol, int model, int interaction, int ngpreds,
//            bool iscox, int nullmodel = 0)
//    {
//        //		cout << "model = " << model << "\n";
//        //		cdata.X.print();
//        coxph_data cdata = cdatain.get_unmasked_data();
//        mematrix<double> X = t_apply_model(cdata.X, model, interaction, ngpreds,
//                iscox, nullmodel);
//        //		X.print();
//        int length_beta = X.nrow;
//        beta.reinit(length_beta, 1);
//        sebeta.reinit(length_beta, 1);
//        mematrix<double> newoffset = cdata.offset;
//        newoffset = cdata.offset - (cdata.offset).column_mean(0);
//        mematrix<double> means(X.nrow, 1);
//        for (int i = 0; i < X.nrow; i++)
//            beta[i] = 0.;
//        mematrix<double> u(X.nrow, 1);
//        mematrix<double> imat(X.nrow, X.nrow);
//        double work[X.ncol * 2 + 2 * (X.nrow) * (X.nrow) + 3 * (X.nrow)];
//        double loglik_int[2];
//        int flag;
//        double sctest = 1.0;
////TODO(maarten): remove the following comment signs. This is done for testing purpose only EVIL!
//        coxfit2(&maxiter, &cdata.nids, &X.nrow, cdata.stime.data,
//                cdata.sstat.data, X.data, newoffset.data, cdata.weights.data,
//                cdata.strata.data, means.data, beta.data, u.data, imat.data,loglik_int, &flag, work, &eps, &tol_chol, &sctest);
//        for (int i = 0; i < X.nrow; i++)
//            sebeta[i] = sqrt(imat.get(i, i));
//        loglik = loglik_int[1];
//        niter = maxiter;
//    }
//};
