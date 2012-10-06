
#include "coxph_reg.h"
#include "coxph_data.h"
extern "C"
{
#include "survproto.h"
}

coxph_reg::coxph_reg(coxph_data &cdatain) {
	coxph_data cdata = cdatain.get_unmasked_data();
	beta.reinit(cdata.X.nrow, 1);
	sebeta.reinit(cdata.X.nrow, 1);
	loglik = -9.999e+32;
	sigma2 = -1.;
	chi2_score = -1.;
	niter = 0;
}
coxph_reg::~coxph_reg() {
	//		delete beta;
	//		delete sebeta;
}
void coxph_reg::estimate(coxph_data &cdatain, int verbose, int maxiter,
		double eps, double tol_chol, int model, int interaction, int ngpreds,
		bool iscox, int nullmodel = 0) {
	//		cout << "model = " << model << "\n";
	//		cdata.X.print();
	coxph_data cdata = cdatain.get_unmasked_data();
	mematrix<double> X = t_apply_model(cdata.X, model, interaction, ngpreds,
			iscox, nullmodel);
	//		X.print();
	int length_beta = X.nrow;
	beta.reinit(length_beta, 1);
	sebeta.reinit(length_beta, 1);
	mematrix<double> newoffset = cdata.offset;
	newoffset = cdata.offset - (cdata.offset).column_mean(0);
	mematrix<double> means(X.nrow, 1);
	for (int i = 0; i < X.nrow; i++)
		beta[i] = 0.;
	mematrix<double> u(X.nrow, 1);
	mematrix<double> imat(X.nrow, X.nrow);
	double work[X.ncol * 2 + 2 * (X.nrow) * (X.nrow) + 3 * (X.nrow)];
	double loglik_int[2];
	int flag;
	double sctest = 1.0;
//TODO(maarten): remove the following comment signs. This is done for testing purpose only EVIL!
	coxfit2(&maxiter, &cdata.nids, &X.nrow, cdata.stime.data, cdata.sstat.data,
			X.data, newoffset.data, cdata.weights.data, cdata.strata.data,
			means.data, beta.data, u.data, imat.data, loglik_int, &flag, work,
			&eps, &tol_chol, &sctest);
	for (int i = 0; i < X.nrow; i++)
		sebeta[i] = sqrt(imat.get(i, i));
	loglik = loglik_int[1];
	niter = maxiter;
}

