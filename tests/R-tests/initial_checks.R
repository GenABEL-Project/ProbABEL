## This file contains the parts of the R code that are shared by the
## three R-based checks for palinear, palogist and pacoxph,
## respectively.

## Set tolerance for comparing various outputs
tol <- 1e-5


####
## load the data
####
example.path <- paste0(srcdir, "../../examples/")

## load phenotypic data
pheno <- read.table(paste0(example.path, pheno.file),
                    head=TRUE, string=FALSE)

## load genetic DOSE data
dose <- read.table(paste0(example.path, "test.mldose"),
                   head=FALSE, string=FALSE)
## remove "1->" from the names of dose-IDs
idNames   <- dose[, 1]
idNames   <- sub("[0-9]+->", "", idNames)
dose[, 1] <- idNames
cat("Dose: check consistency of names\t\t")
stopifnot( all.equal(dose[, 1], pheno[, 1], tol) )
cat("OK\n")

## load genetic PROB data
prob <- read.table(paste0(example.path, "test.mlprob"),
                   head=FALSE, string=FALSE)
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
