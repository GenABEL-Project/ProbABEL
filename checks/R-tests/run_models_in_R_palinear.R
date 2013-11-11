cat("Checking linear regression...\n")

args <- commandArgs(TRUE)
srcdir <- args[1]

if (is.na(srcdir)) {
    srcdir <- "./"
}

pheno.file <- "height.txt"

source(paste0(srcdir, "initial_checks.R"))

####
## Run ProbABEL to get the output data we want to compare/verify
####
cat("Running ProbABEL...\t\t\t\t")
tmp <- system(paste0("cd ", tests.path, "; bash test_qt.sh; cd -"),
              intern=TRUE)
cat("OK\n")

resPaAddDose <- read.table(
    paste0(tests.path, "height_base_add.out.txt"),
    head=TRUE)[, colsAddDose]
resPaAddProb <- read.table(
    paste0(tests.path, "height_ngp2_add.out.txt"),
    head=TRUE)[, colsAddProb]
resPaDom <- read.table(
    paste0(tests.path, "height_ngp2_domin.out.txt"),
    head=TRUE)[, colsDom]
resPaRec <- read.table(
    paste0(tests.path, "height_ngp2_recess.out.txt"),
    head=TRUE)[, colsRec]
resPaOdom <- read.table(
    paste0(tests.path, "height_ngp2_over_domin.out.txt"),
    head=TRUE)[, colsOdom]
resPa2df <- read.table(
    paste0(tests.path, "height_ngp2_2df.out.txt"),
    head=TRUE)[, cols2df]

####
## run analysis in R
####
attach(pheno)

cat("Comparing R output with ProbABEL output\t\t")

source("run_model_linear.R")

model.fn.0 <- "lm( height[noNA] ~ sex[noNA] + age[noNA] )"
model.fn   <- "lm( height ~ sex + age + snp )"

## Additive model, dosages
snpdose <- "dose[, i]"
dose.add.R <- run.model(model.fn.0, model.fn, snpdose)
colnames(dose.add.R) <- colsAddDose
rownames(dose.add.R) <- NULL
stopifnot( all.equal(resPaAddDose, dose.add.R, tol=tol) )
cat("additive ")


## Additive model, probabilities
snpprob <- "doseFromProb[, i]"
prob.add.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.add.R) <- colsAddProb
rownames(prob.add.R) <- NULL
stopifnot( all.equal(resPaAddProb, prob.add.R, tol=tol) )
cat("additive ")

## dominant model
snpprob <- "prob[, indexHom] + prob[, indexHet]"
prob.dom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.dom.R) <- colsDom
rownames(prob.dom.R) <- NULL
stopifnot( all.equal(resPaDom, prob.dom.R, tol=tol) )
cat("dominant ")

## recessive model
snpprob <- "prob[, indexHom]"
prob.rec.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.rec.R) <- colsRec
rownames(prob.rec.R) <- NULL
stopifnot( all.equal(resPaRec, prob.rec.R, tol=tol) )
cat("recessive ")

## over-dominant model
snpprob <- "prob[, indexHet]"
prob.odom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.odom.R) <- colsOdom
rownames(prob.odom.R) <- NULL
stopifnot( all.equal(resPaOdom, prob.odom.R, tol=tol) )
cat("overdominant ")

## 2df model
prob.2df.R <- data.frame()
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- lm( height[noNA] ~ sex[noNA] + age[noNA] )
        model   <- lm( height ~ sex + age + prob[, indexHet] + prob[, indexHom] )
        smA1A2  <- summary(model)$coef[4, 1:2]
        smA1A1  <- summary(model)$coef[5, 1:2]
        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

        rsq <- resPa2df[i-2, "Rsq"]
        if( rsq < rsq.thresh) {
            row <- c(rsq, NaN, NaN, NaN, NaN, NaN)
        } else {
            row <- c(rsq, smA1A2[1], smA1A2[2], smA1A1[1], smA1A1[2], lrt)

        }
        prob.2df.R <- rbind(prob.2df.R, row)
}
colnames(prob.2df.R) <- cols2df
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(resPa2df, prob.2df.R, tol=tol) )
cat("2df\n")

cat("\t\t\t\t\t\tOK\n")
