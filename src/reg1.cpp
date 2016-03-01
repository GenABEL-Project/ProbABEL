/**
 * \file   reg1.cpp
 * \author The ProbABEL team
 *
 * \brief File containing the parts of the code for linear and
 * logistic regression.
 *
 * For CoxPH regression look in the file coxph_data.h.
 */
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


#include "reg1.h"


/**
 * \brief Apply the genetic model to the design matrix X before
 * running the regression.
 *
 * This also includes taking care of an interaction term if the user
 * requested that.
 *
 * Depending on the number of genetic predictors (ngpreds) the integer
 * coding for the models has a different meaning:
 *
 * If ngpreds==1 (dose data):
 * \li model 0 = additive (1 df)
 *
 * If ngpreds==2 (prob data):
 * \li model 0 = 2 df
 * \li model 1 = additive (1 df)
 * \li model 2 = dominant (1 df)
 * \li model 3 = recessive (1 df)
 * \li model 4 = over-dominant (1 df)
 * @param X Design matrix, including SNP column(s).
 * @param model Integer describing the genetic model to be
 * applied. See the list above.
 * @param interaction Column number of the covariate used in the
 * interaction term.
 * @param ngpreds Number of genetic predictors (1 for dosage data, 2
 * for probability data).
 * @param iscox Indicates whether a CoxPH regression is being done.
 * @param nullmodel Indicates whether the null model is being analysed.
 *
 * @return Matrix with the model applied to it.
 */
mematrix<double> apply_model(const mematrix<double>& X,
                             const int model,
                             const int interaction,
                             const int ngpreds,
                             const bool iscox,
                             const int nullmodel)
{
    if (nullmodel)
    {
        // No need to apply any genotypic model when calculating the
        // null model
        return (X);
    }

    if (model == 0) // Only run these steps the first time this
                    // function is called
    {
        if (interaction != 0)
        {
            // The user requested analysis with an interaction term,
            // so ngpreds columns need to be added to the X matrix.
            mematrix<double> nX;
            nX.reinit(X.nrow, X.ncol + ngpreds);
            int csnp_p1 = nX.ncol - 2 * ngpreds;
            int c1 = nX.ncol - ngpreds;
            // The following two variables are only used when ngpreds
            // == 2. Note that the order of the two probabilities is
            // swapped w.r.t. the file (see regdata::update_snp())!
            int csnp_p2 = nX.ncol - 3;
            int c2 = nX.ncol - 1;

            // Copy the data from X to nX (note: nX has more columns!)
            for (int i = 0; i < X.nrow; i++)
                for (int j = 0; j < X.ncol; j++)
                    nX[i * nX.ncol + j] = X[i * X.ncol + j];

            for (int i = 0; i < nX.nrow; i++)
            {
                if (iscox)
                {
                    // Maksim: interaction with SNP;;
                    nX[i * nX.ncol + c1] =
                        X[i * X.ncol + csnp_p1]
                        * X[i * X.ncol + interaction - 1];
                    if (ngpreds == 2)
                    {
                        nX[i * nX.ncol + c2] =
                            X[i * X.ncol + csnp_p2]
                            * X[i * X.ncol + interaction - 1];
                    }
                }
                else
                {
                    // Maksim: interaction with SNP;;
                    nX[i * nX.ncol + c1] =
                        X[i * X.ncol + csnp_p1]
                        * X[i * X.ncol + interaction];
                    if (ngpreds == 2)
                    {
                        nX[i * nX.ncol + c2] =
                            X[i * X.ncol + csnp_p2]
                            * X[i * X.ncol + interaction];
                    }
                }
            }

            return (nX);
        } // End if (interaction !=0)
        else
        {
            // No interaction analysis, no need to add/change columns
            // to/in X.
            return (X);
        }
    } // End: if (model == 0)


    // NOTE: The rest of this function is not run for dosage data
    // (ngpreds == 1) because in that case there is only one value for
    // model: model == 0; so this function stops at the if() that
    // ended just above.
    mematrix<double> nX;
    if (interaction != 0)
    {
        nX.reinit(X.nrow, (X.ncol));
    }
    else
    {
        nX.reinit(X.nrow, (X.ncol - 1));
    }

    // Note that the order of the two probabilities is swapped
    // w.r.t. the file (see regdata::update_snp())!
    // column with Prob(A1A2)
    int c1 = X.ncol - 2;
    // column with Prob(A1A1)
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
        else if (model == 2)  // dominant
            nX[i * nX.ncol + c1] = X[i * X.ncol + c1] + X[i * X.ncol + c2];
        else if (model == 3)  // recessive
            nX[i * nX.ncol + c1] = X[i * X.ncol + c2];
        else if (model == 4)  // over-dominant
            nX[i * nX.ncol + c1] = X[i * X.ncol + c1];

        if (interaction != 0)
            nX[i * nX.ncol + c2] = X[i * nX.ncol + interaction]
                    * nX[i * nX.ncol + c1];  // Maksim: interaction with SNP
    }

    return nX;
}


/**
 * \brief Apply a genetic model to a transposed design matrix
 * \f$X\f$. Similar to apply_model(), but used in case the design
 * matrix is transposed.
 *
 * The function transposes a temporary copy of the input matrix,
 * applies the model (using apply_model()), transposes it back again
 * and returns that matrix.
 *
 * Used only when doing Cox PH regression.
 * @param X The transposed design matrix, including SNP column(s).
 * @param model Integer describing the genetic model to be
 * applied. See apply_model() for details.
 * @param interaction Column number of the covariate used in the
 * interaction term.
 * @param ngpreds Number of genetic predictors (1 for dosage data, 2
 * for probability data).
 * @param iscox Indicates whether a CoxPH regression is being done.
 * @param nullmodel Indicates whether the null model is being analysed.
 *
 * @return (transposed) Matrix with the model applied to it.
 */
mematrix<double> t_apply_model(const mematrix<double>& X,
                               const int model,
                               const int interaction,
                               const int ngpreds,
                               const bool iscox,
                               const int nullmodel)
{
    mematrix<double> tmpX = transpose(X);
    mematrix<double> nX = apply_model(tmpX, model, interaction, ngpreds,
                                      iscox, nullmodel);
    mematrix<double> out = transpose(nX);
    return out;
}


linear_reg::linear_reg(const regdata& rdatain) {
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


void base_reg::base_score(const mematrix<double>& resid,
                          const int model,
                          const int interaction,
                          const int ngpreds,
                          const masked_matrix& invvarmatrix,
                          const int nullmodel) {
    mematrix<double> oX = reg_data.extract_genotypes();
    mematrix<double> X = apply_model(oX, model, interaction, ngpreds,
                                     false, nullmodel);
    beta.reinit(X.ncol, 1);
    sebeta.reinit(X.ncol, 1);
    int length_beta = X.ncol;
    double N = static_cast<double>(resid.nrow);
    mematrix<double> tX = transpose(X);
    if (invvarmatrix.length_of_mask != 0){
        tX = tX * invvarmatrix.masked_data;
    }

    mematrix<double> u = tX * resid;
    mematrix<double> v = tX * X;
    mematrix<double> csum = column_sum(X);
    csum = transpose(csum) * csum;
    csum = csum * (1. / N);
    v = v - csum;
    // use cholesky to invert

    LDLT <MatrixXd> Ch = LDLT < MatrixXd > (v.data.selfadjointView<Lower>());
    // before was
    // mematrix<double> v_i = invert(v);
    beta.data = Ch.solve(v.data.adjoint() * u.data);
    //TODO(maartenk): set size of v_i directly or remove mematrix class
    mematrix<double> v_i = v;
    v_i.data = Ch.solve(MatrixXd(length_beta, length_beta).
                                    Identity(length_beta, length_beta));

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


/**
 * \brief Solve the linear system in case the --mmscore option was
 * specified.
 *
 * Specifying the --mmscore command line option requires a file name
 * argument as well. This file should contain the inverse
 * variance-covariance matrix file. This function is run when Linear
 * regression is done in combination with the mmscore option. It
 * solves the 'mmscore' equation as specified in Eq. (5) in section
 * 8.2.1 of the ProbABEL manual:
 * \f[
 *    \hat{\beta}_g = (\mathbf{X}^T_g
 *    \mathbf{V}^{-1}_{\hat{\sigma}^2,\hat{h}^2}
 *    \mathbf{X}_g)^{-1}
 *    \mathbf{X}^T_g \mathbf{V}^{-1}_{\hat{\sigma}^2,\hat{h}^2}
 *    \mathbf{R}_{\hat{\beta}_x},
 * \f]
 * where \f$\mathbf{V}^{-1}_{\hat{\sigma}^2,\hat{h}^2}\f$ is the
 * inverse variance-covariance matrix, and
 * \f$\mathbf{R}_{\hat{\beta}_x}\f$ is the vector containing the
 * residuals obtained from the base regression model i.e. the
 * phenotype. In this function, the phenotype is stored in the
 * variable \c Y.
 *
 * @param X The design matrix \f$X_g\f$. \c X should only contain the
 * parts involving genotype data (including any interactions involving
 * a genetic term), all other covariates should have been regressed out
 * before running ProbABEL.
 * @param W_masked The inverse variance-covariance matrix
 * \f$\mathbf{V}^{-1}_{\hat{\sigma}^2,\hat{h}^2}\f$.
 * @param Ch Reference to the LDLT Cholesky decomposition of the
 * matrix to be inverted to get \f$\hat\beta_g\f$:
 * \f[
 * \mathbf{X}^T_g
 *    \mathbf{V}^{-1}_{\hat{\sigma}^2,\hat{h}^2}
 *    \mathbf{X}_g.
 * \f]
 * On return this variable contains said matrix.
 */
void linear_reg::mmscore_regression(const mematrix<double>& X,
                                    const masked_matrix& W_masked,
                                    LDLT<MatrixXd>& Ch) {
    VectorXd Y = reg_data.Y.data.col(0);
    /*
     in ProbABEL <0.50 this calculation was performed like t(X)*W
     This changed to W*X since this is better vectorized since the left hand
     side has more rows: this introduces an additional transpose, but can be
     neglected compared to the speedup this brings(about a factor 2 for the
     palinear with 1 predictor)

     This function solves the system
        (X^T W X) beta = X^T W Y.
     Since W is symmetric (WX)^T = X^T W, so this can be rewritten as
        (WX)^T X beta = (WX)^T Y,
     which is solved using LDLT Cholesky decomposition.
     */
    MatrixXd WX = W_masked.masked_data->data * X.data;
    MatrixXd XWT = WX.transpose();
    Ch = LDLT<MatrixXd>(XWT * X.data);
    VectorXd beta_vec = Ch.solve(XWT * Y);
    sigma2 = (Y - WX * beta_vec).squaredNorm();
    beta.data = beta_vec;
}


/**
 * \brief Perform least squares regression
 *
 * Basically we solve the following linear system here:
 * \f[
 * \mathbf{y} = \mathbf{X} \mathbf{\beta}
 * \f]
 *
 * This function also estimates \f$\sigma^2\f$, the variance of the
 * error term. An estimator of \f$\sigma^2\f$ is given by:
 * \f[
 * \hat\sigma^2 = \frac{1}{n-p}||\mathbf{y} - \mathbf{X} \mathbf{\beta} ||^2,
 * \f]
 * with \f$n\f$ the number of rows of \f$\mathbf{X}\f$ and \f$p\f$
 * the number of columns of \f$\mathbf{X}\f$.
 *
 * @param X The design matrix
 * @param Ch
 */
void linear_reg::LeastSquaredRegression(const mematrix<double>& X,
                                        LDLT<MatrixXd>& Ch) {
    int m = X.ncol;
    MatrixXd txx = MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(
            X.data.adjoint());
    Ch = LDLT < MatrixXd > (txx.selfadjointView<Lower>());
    beta.data = Ch.solve(X.data.adjoint() * reg_data.Y.data);
    sigma2 = (reg_data.Y.data - (X.data * beta.data)).squaredNorm();
}


/**
 * \brief Calculate the loglikelihood of a linear regression contained
 * in a linear_reg object.
 *
 * @param X The design matrix.
 */
void linear_reg::logLikelihood(const mematrix<double>& X) {
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
    //loglik -= halfrecsig2 * residuals[i] * residuals[i];

    double intercept = beta.get(0, 0);
    residuals.data = reg_data.Y.data.array() - intercept;
    //matrix.
    ArrayXXd betacol =
            beta.data.block(1, 0, beta.data.rows() - 1, 1).array().transpose();
    ArrayXXd resid_sub = (X.data.block(0, 1, X.data.rows(), X.data.cols() - 1)
            * betacol.matrix().asDiagonal()).rowwise().sum();
    //std::cout << resid_sub << std::endl;
    residuals.data -= resid_sub.matrix();
    //residuals[i] -= resid_sub;
    loglik -= (residuals.data.array().square() * halfrecsig2).sum();
    loglik -= static_cast<double>(reg_data.nids) * log(sqrt(sigma2));
}


/**
 * \brief Calculate the robust standard errors and covariance.
 *
 * @param X The design matrix \f$X\f$.
 * @param robust_sigma2 The variable in which to store the robust
 * \f$\hat{\sigma}^2\f$.
 * @param tXX_inv Matrix containing \f$(X^T X)^{-1}\f$.
 * @param offset Offset indicating where the columns of the SNP
 * parameters start in the linear_reg::beta vector.
 */
void linear_reg::RobustSEandCovariance(const mematrix<double> &X,
                                       mematrix<double> robust_sigma2,
                                       const MatrixXd tXX_inv,
                                       const int offset) {
    MatrixXd Xresiduals = X.data.array().colwise()
        * residuals.data.col(0).array();
    MatrixXd XbyR =
        MatrixXd(X.ncol, X.ncol).setZero().selfadjointView<Lower>().rankUpdate(
            Xresiduals.adjoint());
    robust_sigma2.data = tXX_inv * XbyR * tXX_inv;
    sebeta.data = robust_sigma2.data.diagonal().array().sqrt();
    covariance.data =
        robust_sigma2.data.bottomLeftCorner(offset, offset).diagonal();
}


/**
 * Calculate the (plain) standard error and covariance of the betas
 * of a linear_reg object.
 *
 * @param sigma2_internal The internal \f$\hat{\sigma}^2\f$.
 * @param tXX_inv Matrix containing \f$(X^T X)^{-1}\f$.
 * @param offset Offset indicating where the columns of the SNP
 * parameters start in the linear_reg::beta vector.
 */
void linear_reg::PlainSEandCovariance(const double sigma2_internal,
                                      const MatrixXd &tXX_inv,
                                      const int offset) {
    sebeta.data = (sigma2_internal * tXX_inv.diagonal().array()).sqrt();
    covariance.data = sigma2_internal
            * tXX_inv.bottomLeftCorner(offset, offset).diagonal().array();
}


/**
 * \brief Estimate the parameters for linear regression.
 *
 * @param verbose Turns verbose printing of various matrices on if
 * non-zero.
 * @param model The number of the genetic model (e.g. additive,
 * recessive, ...) that is to be applied by the apply_model() function.
 * @param interaction
 * @param ngpreds Number of genomic predictors (1 for dosages, 2 for
 * probabilities).
 * @param invvarmatrixin
 * @param robust If non-zero calculate robust standard errors.
 * @param nullmodel If non-zero calculate the null model (excluding
 * SNP information).
 */
void linear_reg::estimate(const int verbose,
                          const int model,
                          const int interaction,
                          const int ngpreds,
                          masked_matrix& invvarmatrixin,
                          const int robust,
                          const int nullmodel) {
    // suda interaction parameter
    // model should come here
    //regdata rdata = rdatain.get_unmasked_data();

    if (verbose)
    {
        // std::cout << "invvarmatrix:\n";
        // invvarmatrixin.masked_data->print();
        std::cout << "rdata.X:\n";
        reg_data.X.print();
    }

    mematrix<double> X = apply_model(reg_data.X, model, interaction, ngpreds,
                                     false, nullmodel);
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

    LDLT <MatrixXd> Ch;
    if (invvarmatrixin.length_of_mask != 0)
    {
        //retrieve masked data W
        invvarmatrixin.update_mask(reg_data.masked_data);
        mmscore_regression(X, invvarmatrixin, Ch);
        double N = X.nrow;
        //sigma2_internal = sigma2 / (N - static_cast<double>(length_beta));
        // Ugly fix to the fact that if we do mmscore, sigma2 is already
        //  in the matrix...
        //      YSA, 2009.07.20
        sigma2_internal = 1.0;
        sigma2 /= N;
    }
    else  // NO mm-score regression : normal least square regression
    {
        LeastSquaredRegression(X, Ch);
        double N = static_cast<double>(X.nrow);
        double P = static_cast<double>(length_beta);
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
    logLikelihood(X);

    MatrixXd tXX_inv = Ch.solve(MatrixXd(length_beta, length_beta).
                                Identity(length_beta, length_beta));
    mematrix<double> robust_sigma2(X.ncol, X.ncol);

    int offset = X.ncol- 1;
     //if additive and interaction and 2 predictors and more than 2 betas
     if (model == 0 && interaction != 0 && ngpreds == 2 && length_beta > 2){
         offset = X.ncol - 2;
     }

    if (robust)
    {
        RobustSEandCovariance(X, robust_sigma2, tXX_inv, offset);
    }
    else
    {
        PlainSEandCovariance(sigma2_internal, tXX_inv, offset);
    }
}


void linear_reg::score(const mematrix<double>& resid,
                       const int model,
                       const int interaction, const int ngpreds,
                       const masked_matrix& invvarmatrix,
                       int nullmodel)
{
    //regdata rdata = rdatain.get_unmasked_data();
    base_score(resid, model, interaction, ngpreds,
               invvarmatrix, nullmodel = 0);
    // TODO: find out why nullmodel is assigned 0 in the call above.
}


logistic_reg::logistic_reg(const regdata& rdatain)
{
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


void logistic_reg::estimate(const int verbose,
                            const int model,
                            const int interaction, const int ngpreds,
                            masked_matrix& invvarmatrixin,
                            const int robust, const int nullmodel) {
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
                                     false, nullmodel);
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
        //TODO(maarten) invvarmatix is symmetric:is there an more effective way?
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
    while (niter < MAXITER && delta > EPS)
    {
        mematrix<double> eMu = (X) * beta;
        mematrix<double> eMu_us = eMu;
        for (int i = 0; i < eMu.nrow; i++)
        {
            double emu = eMu.get(i, 0);
            double value = emu;
            double zval;
            double expval = exp(value);
            value = expval / (1. + expval);
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

        // Compute the likelihood.
        double prevlik = loglik;
        loglik = 0;
        for (int i = 0; i < eMu.nrow; i++) {
            loglik += reg_data.Y[i] * eMu_us[i] - log(1. + exp(eMu_us[i]));
        }

        delta = fabs(1. - (prevlik / loglik));
        niter++;
    }  // END while (niter < MAXITER && delta > EPS)

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


void logistic_reg::score(const mematrix<double>& resid,
                         const int model,
                         const int interaction, const int ngpreds,
                         const masked_matrix& invvarmatrix,
                         int nullmodel) {
    base_score(resid, model, interaction, ngpreds,
               invvarmatrix, nullmodel = 0);
    // TODO: find out why nullmodel is assigned 0 in the call above.
}
