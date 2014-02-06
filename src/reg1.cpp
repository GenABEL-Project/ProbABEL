#include "reg1.h"

mematrix<double> apply_model(mematrix<double>& X, int model, int interaction,
        int ngpreds, bool is_interaction_excluded, bool iscox, int nullmodel)
// if ngpreds==1 (dose data):
// model 0 = additive 1 df
// if ngpreds==2 (prob data):
// model 0 = 2 df
// model 1 = additive 1 df
// model 2 = dominant 1 df
// model 3 = recessive 1 df
// model 4 = over-dominant 1 df
        {
    if (nullmodel)
    {
        // No need to apply any genotypic model when calculating the
        // null model
        return (X);
    }

    if (model == 0)
    {
        if (interaction != 0 && !nullmodel)
        {
            if (ngpreds == 2)
            {
                mematrix<double> nX;
                nX.reinit(X.nrow, X.ncol + 2);
                int csnp_p1 = nX.ncol - 4;
                int csnp_p2 = nX.ncol - 3;
                int c1 = nX.ncol - 2;
                int c2 = nX.ncol - 1;
                for (int i = 0; i < X.nrow; i++)
                    for (int j = 0; j < (X.ncol); j++)
                        nX[i * nX.ncol + j] = X[i * X.ncol + j];
                for (int i = 0; i < nX.nrow; i++)
                {
                    if (iscox)
                    {
                        // Maksim: interaction with SNP;;
                        nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
                                * X[i * X.ncol + interaction - 1];
                        nX[i * nX.ncol + c2] = X[i * X.ncol + csnp_p2]
                                * X[i * X.ncol + interaction - 1];
                    }
                    else
                    {
                        // Maksim: interaction with SNP;;
                        nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
                                * X[i * X.ncol + interaction];
                        nX[i * nX.ncol + c2] = X[i * X.ncol + csnp_p2]
                                * X[i * X.ncol + interaction];
                    }
                }
                //________________________
                if (is_interaction_excluded)
                {
                    mematrix<double> nX_without_interact_phe;
                    nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);

                    for (int row = 0; row < nX.nrow; row++)
                    {
                        // Han Chen
                        int col_new = -1;
                        for (int col = 0; col < nX.ncol; col++)
                        {
                            if (col != interaction && !iscox)
                            {
                                col_new++;
                                nX_without_interact_phe[row
                                        * nX_without_interact_phe.ncol + col_new] =
                                        nX[row * nX.ncol + col];
                            }
                            if (col != interaction - 1 && iscox)
                            {
                                col_new++;
                                nX_without_interact_phe[row
                                        * nX_without_interact_phe.ncol + col_new] =
                                        nX[row * nX.ncol + col];
                            }
                        } // interaction_only, model==0, ngpreds==2
                          // Oct 26, 2009
                    }
                    return nX_without_interact_phe;
                }  // end of is_interaction_excluded
                   //________________________
                return (nX);
            }
            if (ngpreds == 1)
            {
                mematrix<double> nX;
                nX.reinit(X.nrow, X.ncol + 1);
                int csnp_p1 = nX.ncol - 2;
                int c1 = nX.ncol - 1;
                for (int i = 0; i < X.nrow; i++)
                    for (int j = 0; j < (X.ncol); j++)
                        nX[i * nX.ncol + j] = X[i * X.ncol + j];

                for (int i = 0; i < nX.nrow; i++)
                {
                    if (iscox)
                    {
                        // Maksim: interaction with SNP;;
                        nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
                                * X[i * X.ncol + interaction - 1];
                    }
                    else
                    {
                        // Maksim: interaction with SNP;;
                        nX[i * nX.ncol + c1] = X[i * X.ncol + csnp_p1]
                                * X[i * X.ncol + interaction];
                    }
                }
                //________________________
                if (is_interaction_excluded)
                {
                    mematrix<double> nX_without_interact_phe;
                    nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);
                    for (int row = 0; row < nX.nrow; row++)
                    {
                        int col_new = -1;
                        for (int col = 0; col < nX.ncol; col++)
                        {
                            if (col != interaction && !iscox)
                            {
                                col_new++;
                                nX_without_interact_phe[row
                                        * nX_without_interact_phe.ncol + col_new] =
                                        nX[row * nX.ncol + col];
                            }
                            if (col != interaction - 1 && iscox)
                            {
                                col_new++;
                                nX_without_interact_phe[row
                                        * nX_without_interact_phe.ncol + col_new] =
                                        nX[row * nX.ncol + col];
                            }
                        }
                    }
                    return nX_without_interact_phe;
                }  // end of is_interaction_excluded
                   //________________________
                return (nX);
            }
        }
        else
        {
            return (X);
        }
    }

    mematrix<double> nX;
    if (interaction != 0)
    {
        nX.reinit(X.nrow, (X.ncol));
    }
    else
    {
        nX.reinit(X.nrow, (X.ncol - 1));
    }

    // column with Prob(A1A2)
    int c1 = X.ncol - 2;
    // column with Prob(A1A1). Note the order is swapped cf the file!
    int c2 = X.ncol - 1;

    for (int i = 0; i < X.nrow; i++){
        for (int j = 0; j < (X.ncol - 2); j++){
            nX[i * nX.ncol + j] = X[i * X.ncol + j];
        }
    }

    for (int i = 0; i < nX.nrow; i++)
    {
        if (model == 1)  // additive
            nX[i * nX.ncol + c1] = X[i * X.ncol + c1] + 2. * X[i * X.ncol + c2];
        else if (model == 2)  //dominant
            nX[i * nX.ncol + c1] = X[i * X.ncol + c1] + X[i * X.ncol + c2];
        else if (model == 3)  // recessive
            nX[i * nX.ncol + c1] = X[i * X.ncol + c2];
        else if (model == 4)  // over-dominant
            nX[i * nX.ncol + c1] = X[i * X.ncol + c1];

        if (interaction != 0)
            nX[i * nX.ncol + c2] = X[i * nX.ncol + interaction]
                    * nX[i * nX.ncol + c1];  // Maksim: interaction with SNP
    }

    //Han Chen
    if (is_interaction_excluded)
    {
        mematrix<double> nX_without_interact_phe;
        nX_without_interact_phe.reinit(nX.nrow, nX.ncol - 1);

        for (int row = 0; row < nX.nrow; row++)
        {
            int col_new = -1;
            for (int col = 0; col < nX.ncol; col++)
            {
                if (col != interaction && !iscox)
                {
                    col_new++;
                    nX_without_interact_phe[row * nX_without_interact_phe.ncol
                            + col_new] = nX[row * nX.ncol + col];
                }
                if (col != interaction - 1 && iscox)
                {
                    col_new++;
                    nX_without_interact_phe[row * nX_without_interact_phe.ncol
                            + col_new] = nX[row * nX.ncol + col];
                }
            }
        }
        return nX_without_interact_phe;
    }  // interaction_only, model!=0, ngpreds==2
    return nX;
}

mematrix<double> t_apply_model(mematrix<double>& X, int model, int interaction,
        int ngpreds, bool iscox, int nullmodel) {
    mematrix<double> tmpX = transpose(X);
    mematrix<double> nX = apply_model(tmpX, model, interaction, ngpreds,
            interaction, iscox, nullmodel);
    mematrix<double> out = transpose(nX);
    return out;
}

linear_reg::linear_reg(regdata& rdatain) {
    reg_data = rdatain.get_unmasked_data();
    // std::cout << "linear_reg: " << rdata.nids << " " << (rdata.X).ncol
    //           << " " << (rdata.Y).ncol << "\n";
    int length_beta = (reg_data.X).ncol;
    beta.reinit(length_beta, 1);
    sebeta.reinit(length_beta, 1);
    //Han Chen
    if (length_beta > 1)
    {
        covariance.reinit(length_beta - 1, 1);
    }
    //Oct 26, 2009
    residuals.reinit(reg_data.nids, 1);
    sigma2 = -1.;
    loglik = -9.999e+32;
    chi2_score = -1.;
}

void base_reg::base_score(mematrix<double>& resid,
        double tol_chol, int model, int interaction, int ngpreds,
        const masked_matrix& invvarmatrix, int nullmodel) {
    mematrix<double> oX = reg_data.extract_genotypes();
    mematrix<double> X = apply_model(oX, model, interaction, ngpreds,
            reg_data.is_interaction_excluded, false, nullmodel);
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
    for (int i = 0; i < resid.nrow; i++)
    {
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


void linear_reg::estimate( int verbose, double tol_chol,
        int model, int interaction, int ngpreds, masked_matrix& invvarmatrixin,
        int robust, int nullmodel) {
    // suda interaction parameter
    // model should come here
    //regdata rdata = rdatain.get_unmasked_data();

    if (verbose)
    {
        cout << reg_data.is_interaction_excluded
                << " <-rdata.is_interaction_excluded\n";
        // std::cout << "invvarmatrix:\n";
        // invvarmatrixin.masked_data->print();
        std::cout << "rdata.X:\n";
        reg_data.X.print();
    }

    mematrix<double> X = apply_model(reg_data.X, model, interaction, ngpreds,
            reg_data.is_interaction_excluded, false, nullmodel);
    if (verbose)
    {
        std::cout << "X:\n";
        X.print();
        std::cout << "Y:\n";
        reg_data.Y.print();
    }

    int length_beta = X.ncol;
    beta.reinit(length_beta, 1);
    sebeta.reinit(length_beta, 1);
    //Han Chen
    if (length_beta > 1)
    {
        if (model == 0 && interaction != 0 && ngpreds == 2 && length_beta > 2)
        {
            covariance.reinit(length_beta - 2, 1);
        }
        else
        {
            covariance.reinit(length_beta - 1, 1);
        }
    }

    double sigma2_internal;
    mematrix<double> tXX_i;
#if EIGEN
    LDLT <MatrixXd> Ch ;
#endif
    if (invvarmatrixin.length_of_mask != 0)
    {
        //retrieve masked data W
        invvarmatrixin.update_mask(reg_data.masked_data);

        // This regression is Weighted Least Square: used for mmscore :
        // FLOPS count are calculated for 3*1000 matrix as follow:
        //C=AB (m X n matrix A and n x P matrix B)
        //flops=mp(2n-1) (when n is big enough flops=mpn2)
        //Oct 26, 2009
        mematrix<double> tXW = transpose(X) * invvarmatrixin.masked_data; // flops 5997000
        tXX_i = tXW * X;        // 17991 flops
#if EIGEN
        Ch=LDLT <MatrixXd>(tXX_i.data.selfadjointView<Lower>());
#endif

        // use cholesky to invert
        cholesky2_mm(tXX_i, tol_chol);
        chinv2_mm(tXX_i);
        beta = tXX_i * (tXW * reg_data.Y);        // flops 15+5997
        // now compute residual variance
        sigma2 = 0.;
        mematrix<double> sigma2_matrix = reg_data.Y - (transpose(tXW) * beta); //flops: 1000+5000
        for (int i = 0; i < sigma2_matrix.nrow; i++)
        {
            double val = sigma2_matrix.get(i, 0);
            sigma2 += val * val; // flops: 3000
        }
        double N = X.nrow;
        //sigma2_internal = sigma2 / (N - static_cast<double>(length_beta));
        // Ugly fix to the fact that if we do mmscore, sigma2 is already
        //  in the matrix...
        //      YSA, 2009.07.20
        sigma2_internal = 1.0;
        sigma2 /= N;

    }
    else//NO mm-score regression : normal least square regression
    {
#if EIGEN
        int m = X.ncol;
        MatrixXd txx = MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(X.data.adjoint());
        Ch=LDLT <MatrixXd>(txx.selfadjointView<Lower>());
        beta.data= Ch.solve(X.data.adjoint() * reg_data.Y.data);

       tXX_i.data=Ch.solve(MatrixXd(m, m).Identity(m,m));
       tXX_i.nrow=tXX_i.data.rows();
       tXX_i.ncol=tXX_i.data.cols();
       tXX_i.nelements=tXX_i.ncol*tXX_i.nrow;

#else
        mematrix<double> tX = transpose(X);
        // use cholesky to invert
                tXX_i = tX * X;
                cholesky2_mm(tXX_i, tol_chol);
                chinv2_mm(tXX_i);
                beta = tXX_i * (tX * (reg_data.Y));
#endif

        // now compute residual variance
        sigma2 = 0.;
        mematrix<double> sigma2_matrix = reg_data.Y - (X * beta);
#if EIGEN
        sigma2 = sigma2_matrix.data.squaredNorm() ;
#else
        for (int i = 0; i < sigma2_matrix.nrow; i++)
        {
            double val = sigma2_matrix.get(i, 0);
            sigma2 += val * val;
        }
#endif
        double N = static_cast<double>(X.nrow);
        double P=static_cast<double>(length_beta);
        sigma2_internal = sigma2 / (N - P);
        sigma2 /= N;
    }
    /*
     loglik = 0.;
     double ss=0;
     for (int i=0;i<rdata.nids;i++) {
     double resid = rdata.Y[i] - beta.get(0,0); // intercept
     for (int j=1;j<beta.nrow;j++) resid -= beta.get(j,0)*X.get(i,j);
     // residuals[i] = resid;
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


#if EIGEN
    double intercept = beta.get(0, 0);
    residuals.data= reg_data.Y.data.array()-intercept;
    //matrix.
    ArrayXXd betacol = beta.data.block(1,0,beta.data.rows()-1,1).array().transpose();
    ArrayXXd resid_sub = (X.data.block(0,1,X.data.rows(),X.data.cols()-1)*betacol.matrix().asDiagonal()).rowwise().sum() ;
    //std::cout << resid_sub << std::endl;
    residuals.data-=resid_sub.matrix();
    //residuals[i] -= resid_sub;
    loglik-=(residuals.data.array().square()*halfrecsig2).sum();

    //loglik -= halfrecsig2 * residuals[i] * residuals[i];

#else
    for (int i = 0; i < reg_data.nids; i++)
     {
         double resid = reg_data.Y[i] - beta.get(0, 0); // intercept
         for (int j = 1; j < beta.nrow; j++){
             resid -= beta.get(j, 0) * X.get(i, j);
         }
         residuals[i] = resid;
         loglik -= halfrecsig2 * resid * resid;
     }
#endif

    loglik -= static_cast<double>(reg_data.nids) * log(sqrt(sigma2));
#if EIGEN
    MatrixXd tXX_inv=Ch.solve(MatrixXd(length_beta, length_beta).Identity(length_beta,length_beta));
#endif

    mematrix<double> robust_sigma2(X.ncol, X.ncol);
    if (robust)
    {
#if EIGEN
        MatrixXd Xresiduals = X.data.array().colwise()*residuals.data.col(0).array();
        MatrixXd  XbyR = MatrixXd(X.ncol, X.ncol).setZero().selfadjointView<Lower>().rankUpdate(Xresiduals.adjoint());
        robust_sigma2.data= tXX_inv*XbyR *tXX_inv;
#else

        mematrix<double> XbyR = X;
        for (int i = 0; i < X.nrow; i++){
            for (int j = 0; j < X.ncol; j++)
            {
                double tmpval = XbyR.get(i, j) * residuals[i];
                XbyR.put(tmpval, i, j);
            }
        }
        XbyR = transpose(XbyR) * XbyR;
        robust_sigma2 = tXX_i * XbyR;
        robust_sigma2 = robust_sigma2 * tXX_i;

#endif


    }
    //cout << "estimate 0\n";
#if EIGEN
    if (robust)
    {
        sebeta.data = robust_sigma2.data.diagonal().array().sqrt();
    }
    else
    {
        sebeta.data =
                (sigma2_internal
                        * tXX_inv.diagonal().array()).sqrt();
    }
    int offset=X.ncol- 1;
    //if additive and interaction and 2 predictors and more then 2 betas

    if (model == 0 && interaction != 0 && ngpreds == 2 && length_beta > 2){
          offset=X.ncol - 2;
    }

    if (robust)
    {
        covariance.data = robust_sigma2.data.bottomLeftCorner(
                offset, offset).diagonal();

    }
    else
    {
            covariance.data = sigma2_internal
                    * tXX_inv.bottomLeftCorner(offset,
                            offset).diagonal().array();
        }

#else

    //cout << "estimate 0\n";
    for (int i = 0; i < (length_beta); i++)
    {
        if (robust)
        {
            // cout << "estimate :robust\n";
            double value = sqrt(robust_sigma2.get(i, i));
            sebeta.put(value, i, 0);
            //Han Chen
            if (i > 0)
            {
                if (model == 0 && interaction != 0 && ngpreds == 2
                        && length_beta > 2)
                {
                    if (i > 1)
                    {
                        double covval = robust_sigma2.get(i, i - 2);
                        covariance.put(covval, i - 2, 0);
                    }
                }
                else
                {
                    double covval = robust_sigma2.get(i, i - 1);
                    covariance.put(covval, i - 1, 0);
                }
            }
            //Oct 26, 2009
        }
        else
        {
            // cout << "estimate :non-robust\n";
            double value = sqrt(sigma2_internal * tXX_i.get(i, i));
            sebeta.put(value, i, 0);
            //Han Chen
            if (i > 0)
            {
                if (model == 0 && interaction != 0 && ngpreds == 2
                        && length_beta > 2)
                {
                    if (i > 1)
                    {
                        double covval = sigma2_internal * tXX_i.get(i, i - 2);
                        covariance.put(covval, i - 2, 0);
                    }
                }
                else
                {
                    double covval = sigma2_internal * tXX_i.get(i, i - 1);
                    covariance.put(covval, i - 1, 0);
                }
            }
            //Oct 26, 2009
        }
    }
#endif

}

void linear_reg::score(mematrix<double>& resid,
        double tol_chol, int model, int interaction, int ngpreds,
        const masked_matrix& invvarmatrix, int nullmodel) {
   // regdata rdata = rdatain.get_unmasked_data();
    base_score(resid,  tol_chol, model, interaction, ngpreds,
            invvarmatrix, nullmodel = 0);
}

logistic_reg::logistic_reg(regdata& rdatain) {
    reg_data = rdatain.get_unmasked_data();
    int length_beta = reg_data.X.ncol;
    beta.reinit(length_beta, 1);
    sebeta.reinit(length_beta, 1);
    //Han Chen
    if (length_beta > 1)
    {
        covariance.reinit(length_beta - 1, 1);
    }
    //Oct 26, 2009
    residuals.reinit(reg_data.X.nrow, 1);
    sigma2 = -1.;
    loglik = -9.999e+32; // should actually be MAX of the corresponding type
    niter = -1;
    chi2_score = -1.;
}

void logistic_reg::estimate( int verbose, int maxiter,
        double eps, int model, int interaction, int ngpreds,
        masked_matrix& invvarmatrixin, int robust, int nullmodel) {
    // In contrast to the 'linear' case 'invvarmatrix' contains the
    // inverse of correlation matrix (not the inverse of var-cov matrix)
    // h2.object$InvSigma * h.object2$h2an$estimate[length(h2$h2an$estimate)]
    // the inverse of var-cov matrix scaled by total variance
    //regdata rdata = rdatain.get_unmasked_data();
    // a lot of code duplicated between linear and logistic...
    // e.g. a piece below...
    mematrix<double> invvarmatrix;
    if (invvarmatrixin.length_of_mask != 0)
    {
        invvarmatrixin.update_mask(reg_data.masked_data);
    }
    mematrix<double> X = apply_model(reg_data.X, model, interaction, ngpreds,
            reg_data.is_interaction_excluded, false, nullmodel);
    int length_beta = X.ncol;
    beta.reinit(length_beta, 1);
    sebeta.reinit(length_beta, 1);
    //Han Chen

    if (length_beta > 1)
    {

        if (model == 0 && interaction != 0 && ngpreds == 2 && length_beta > 2)
        {
            covariance.reinit(length_beta - 2, 1);
        }
        else
        {
            covariance.reinit(length_beta - 1, 1);
        }
    }

    //Oct 26, 2009
    mematrix<double> W((X).nrow, 1);
    mematrix<double> z((X).nrow, 1);
    mematrix<double> tXWX(length_beta, length_beta);
    mematrix<double> tXWX_i(length_beta, length_beta);
    mematrix<double> tXWz(length_beta, 1);
    double prev = (reg_data.Y).column_mean(0);
    if (prev >= 1. || prev <= 0.)
    {
        std::cerr << "prevalence not within (0,1)\n";
        exit(1);
    }
    for (int i = 0; i < length_beta; i++)
        beta.put(0., i, 0);

    beta.put(log(prev / (1. - prev)), 0, 0);
    mematrix<double> tX = transpose(X);

    if (invvarmatrix.nrow != 0 && invvarmatrix.ncol != 0)
    {
        //TODO(maarten):invvarmatix is symmetric:is there an more effective way?
        tX = tX * invvarmatrix;
    }
    /*
     std::cout << "\n";
     std::cout << "X " << X.get(0,0) << " " << X.get(0,1) << " "
     << X.get(0,2) << "\n";
     if (X.ncol==4) std::cout << "X[4] " << X.get(0,3) << "\n";
     std::cout << "Inv " << invvarmatrix.get(0,0) << " "
     << invvarmatrix.get(0,1) << " "
     << invvarmatrix.get(0,2) << "\n";

     if (X.ncol==4) std::cout << ,"X[4] " << invvarmatrix.get(0,3) << "\n";
     std::cout << "tXInv " << tX.get(0,0) << " " << tX.get(1,0) << " "
     << tX.get(2,0) << "%f\n";
     if (X.ncol==4) std::cout << "X[4] " << tX.get(3,0) << "\n";
     */
    niter = 0;
    double delta = 1.;
    double prevlik = 0.;
    while (niter < maxiter && delta > eps)
    {
        mematrix<double> eMu = (X) * beta;
        mematrix<double> eMu_us = eMu;
        for (int i = 0; i < eMu.nrow; i++)
        {
            double emu = eMu.get(i, 0);
            double value = emu;
            double zval;
            value = exp(value) / (1. + exp(value));
            residuals[i] = (reg_data.Y).get(i, 0) - value;
            eMu.put(value, i, 0);
            W.put(value * (1. - value), i, 0);
            zval = emu
                    + (1. / (value * (1. - value)))
                            * (((reg_data.Y).get(i, 0)) - value);
            z.put(zval, i, 0);
        }

        mematrix<double> tmp = productMatrDiag(tX, W);
        if (verbose)
        {
            std::cout << "tXW:\n";
            tmp.print();
        }
        mematrix<double> tXWX = tmp * (X);
        //N = tXWX.get(0, 0);
        if (verbose)
        {
            std::cout << "tXWX:\n";
            tXWX.print();
        }
        // std::cout << "tXWX:\n";tXWX.print();
        //
        // use cholesky to invert
        //
        // tXWX_i = tXWX;
        //cholesky2_mm(tXWX_i,tol_chol);
        //if (verbose) {std::cout << "chole tXWX:\n"; tXWX_i.print();}
        //std::cout << "chole tXWX:\n"; tXWX_i.print();
        //chinv2_mm(tXWX_i);
        // was before
        tXWX_i = invert(tXWX);
        if (verbose)
        {
            std::cout << "tXWX-1:\n";
            tXWX_i.print();
        }
        // std::cout << "*** tXWX_i\n"; tXWX_i.print();
        mematrix<double> tmp1 = productMatrDiag(tX, W);
        mematrix<double> tXWz = tmp1 * z;
        if (verbose)
        {
            std::cout << "tXWz:\n";
            tXWz.print();
        }
        beta = tXWX_i * tXWz;
        // std::cout << "*** res: " << residuals[0] << " "
        //           << residuals[1] << " " << residuals[2] << "\n";
        //mematrix<double> txres = tx * residuals;
        // std::cout << "*** txres\n";txres.print();
        //beta = txwx_i* txres;
        if (verbose)
        {
            std::cout << "beta:\n";
            beta.print();
        }
        // std::cout << "beta:\n"; beta.print();
        // compute likelihood
        prevlik = loglik;
        loglik = 0.;
        for (int i = 0; i < eMu.nrow; i++)
            loglik += reg_data.Y[i] * eMu_us[i] - log(1. + exp(eMu_us[i]));

        delta = fabs(1. - (prevlik / loglik));
        niter++;
    }  // END while (niter < maxiter && delta > eps)

    sigma2 = 0.;
    mematrix<double> robust_sigma2(X.ncol, X.ncol);

    if (robust)
    {
        mematrix<double> XbyR = X;
        for (int i = 0; i < X.nrow; i++)
            for (int j = 0; j < X.ncol; j++)
            {
                double tmpval = XbyR.get(i, j) * residuals[i];
                XbyR.put(tmpval, i, j);
            }
        XbyR = transpose(XbyR) * XbyR;
        robust_sigma2 = tXWX_i * XbyR;
        robust_sigma2 = robust_sigma2 * tXWX_i;
    }
    for (int i = 0; i < (length_beta); i++)
    {
        if (robust)
        {
            double value = sqrt(robust_sigma2.get(i, i));
            sebeta.put(value, i, 0);
            //Han Chen
            if (i > 0)
            {
                if (model == 0 && interaction != 0 && ngpreds == 2
                        && length_beta > 2)
                {
                    if (i > 1)
                    {
                        double covval = robust_sigma2.get(i, i - 2);
                        covariance.put(covval, i - 2, 0);
                    }
                }
                else
                {
                    double covval = robust_sigma2.get(i, i - 1);
                    covariance.put(covval, i - 1, 0);
                }
            }
            //Oct 26, 2009
        }
        else
        {
            double value = sqrt(tXWX_i.get(i, i));
            sebeta.put(value, i, 0);
            //Han Chen
            if (i > 0)
            {
                if (model == 0 && interaction != 0 && ngpreds == 2
                        && length_beta > 2)
                {
                    if (i > 1)
                    {
                        double covval = tXWX_i.get(i, i - 2);
                        covariance.put(covval, i - 2, 0);
                    }
                }
                else
                {
                    double covval = tXWX_i.get(i, i - 1);
                    covariance.put(covval, i - 1, 0);
                }
            }
            //Oct 26, 2009
        }
    }
    if (verbose)
    {
        std::cout << "sebeta (" << sebeta.nrow << "):\n";
        sebeta.print();
    }
    // std::cout << "beta (" << beta.nrow << "):\n"; beta.print();
    // std::cout << "sebeta (" << sebeta.nrow << "):\n"; sebeta.print();
    // exit(1);
}

void logistic_reg::score(mematrix<double>& resid,
        double tol_chol, int model, int interaction, int ngpreds,
        masked_matrix& invvarmatrix, int nullmodel) {
    base_score(resid, tol_chol, model, interaction, ngpreds,
            invvarmatrix, nullmodel = 0);
}
