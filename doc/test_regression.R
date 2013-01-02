
library(MASS)
NCOV = 2
NOBS <- 1000
beta <- c(0,0,.1)
rvar <- 1
X <- matrix(c(rep(1,NOBS),rnorm(NCOV*NOBS)),ncol=NCOV+1)
Y <- X %*% beta + rnorm(NOBS,sd=sqrt(rvar))
lmf <- lm(Y ~ X[,c(2:(NCOV+1))])
slmf <- summary(lmf)

slmf$sigma
slmf$cov
sqrt(slmf$sigma*diag(slmf$cov)*NOBS/(NOBS-NCOV-1))

XpX <- t(X) %*% X
# == slmf$cov
XpX_i <- ginv(XpX)
slmf$cov
XpX_i
# estimate of beta == smlf$coef
estbeta <- XpX_i %*% (t(X) %*% Y)
slmf$coef
estbeta
# residuals == slmf$resid
resid <- Y - X %*% estbeta
slmf$resid[1:10]
resid[1:10]
# residual variance == slmf$sigma^2
estsigma2 <- as.numeric(t(resid) %*% resid / (NOBS-NCOV-1))
slmf$sigma^2
estsigma2
# inverse of var/covar matrix for parame estimates
estvc <- estsigma2 * XpX_i
slmf$sigma*slmf$cov
estvc
# standard errors
estse <- sqrt(diag(estvc))
slmf$coef[,2]
slmf$sigma*sqrt(diag(slmf$cov))
estse

# Wald test on covariates
# inverse estvc
estvc_i <- ginv(estvc)
estvc_i
1/estse
1/sqrt(diag(estvc))
sqrt(diag(estvc_i))
sqrt(diag(ginv(estvc[c(2:(NCOV+1)),c(2:(NCOV+1))])))

testbeta <- estbeta
testbeta[2:(NCOV+1)] <- 0
t2 <- t(estbeta-testbeta) %*% estvc_i %*% (estbeta-testbeta)
t2
pchisq(t2,NCOV,low=F)


null <- lm(Y~1)
anova(null,lmf,test="Chisq")

fisherI <- estvc_i[c(2,3,1),c(2,3,1)]
I11 <- fisherI[1:2,1:2]
I12 <- fisherI[1:2,3]
I21 <- fisherI[3,1:2]
I22 <- fisherI[3,3]
I1 <- I11 - I12 %*% ginv(I22) %*% I21
t2 <- t(estbeta[2:(NCOV+1)]) %*% I1 %*% estbeta[2:(NCOV+1)]
t2
pchisq(t2,NCOV,low=F)

