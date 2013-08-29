## This file contains the parts of the R code that are shared by the
## three R-based checks for palinear, palogist and pacoxph,
## respectively.

## Set tolerance for comparing various outputs
tol <- 1e-5

## R^2 threshold from the mlinfo file. If R^2 smaller than this
## threshold, beta, se_beta and chi^2 are set to NaN.
## Check this value against the one in the ProbABEL code (main.cpp),
## look for the variables freq and poly.
rsq.thresh <- 1e-16

####
## load the data
####
inputfiles.path <- paste0(srcdir, "../inputfiles/")
tests.path      <- paste0(srcdir, "../")

## load phenotypic data
pheno <- read.table(paste0(inputfiles.path, pheno.file),
                    head=TRUE, string=FALSE)

## load genetic DOSE data
dose <- read.table(paste0(inputfiles.path, "test.mldose"),
                   head=FALSE, string=FALSE)
## remove "1->" from the names of dose-IDs
idNames   <- dose[, 1]
idNames   <- sub("[0-9]+->", "", idNames)
dose[, 1] <- idNames
cat("Dose: check consistency of names\t\t")
stopifnot( all.equal(dose[, 1], pheno[, 1], tol) )
cat("OK\n")

## load genetic PROB data
prob <- read.table(paste0(inputfiles.path, "test.mlprob"),
                   header=FALSE, stringsAsFactors=FALSE)
## remove "1->" from the names of prob-IDs
idNames   <- prob[, 1]
idNames   <- sub("[0-9]+->", "", idNames)
prob[, 1] <- idNames
cat("Prob: check consistency of names\t\t")
stopifnot( all.equal(prob[, 1], pheno[, 1], tol) )
cat("OK\n")

## check consistency DOSE <-> PROB
doseFromProb <- matrix(NA, ncol=dim(dose)[2], nrow=dim(dose)[1])
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        doseFromProb[, i] <- prob[, indexHom] * 2 + prob[, indexHet]
}
cat("Check consistency dose <-> prob gtdata\t\t")
stopifnot( all.equal(dose[, 3:ncol(dose)],
                     as.data.frame(doseFromProb)[,3:ncol(doseFromProb)],
                     tol=tol )
          )
cat("OK\n")


## Read the imputed Rsq from the info file
Rsq <- read.table(paste0(inputfiles.path, "test.mlinfo"),
                  header=TRUE,
                  stringsAsFactors=FALSE)[, c("Rsq")]


## Define column names of the various ProbABEL output file headers
colsAddDose <- c("Rsq",
                 "beta_SNP_add",
                 "sebeta_SNP_add",
                 "chi2_SNP")
colsAddProb <- c("Rsq",
                 "beta_SNP_addA1",
                 "sebeta_SNP_addA1",
                 "chi2_SNP_A1")
colsDom <- c("Rsq",
             "beta_SNP_domA1",
             "sebeta_SNP_domA1",
             "chi2_SNP_domA1")
colsRec <- c("Rsq",
             "beta_SNP_recA1",
             "sebeta_SNP_recA1",
             "chi2_SNP_recA1")
colsOdom <-c("Rsq",
             "beta_SNP_odomA1",
             "sebeta_SNP_odomA1",
             "chi2_SNP_odomA1")
cols2df <- c("Rsq",
             "beta_SNP_A1A2",
             "sebeta_SNP_A1A2",
             "beta_SNP_A1A1",
             "sebeta_SNP_A1A1",
             "chi2_SNP_2df")
