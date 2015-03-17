cat("Checking Cox PH regression...\n")
if (!require(survival)) {
    cat("The R package 'survival' is not installed. Skipping Cox PH checks\n")
    q()
}

args <- commandArgs(TRUE)
srcdir <- args[1]

if (is.na(srcdir)) {
    srcdir <- "./"
} else {
    ## Apparently we are running R from the command line. Disable
    ## warnings so that they don't clutter the screen when running
    ## this script.
    old.warn <- options()$warn
    options(warn=-1)
}


pheno.file <- "coxph_data.txt"

source(paste0(srcdir, "initial_checks.R"))

####
## Run ProbABEL to get the output data we want to compare/verify
####
cat("Running ProbABEL...\t\t\t\t")
tmp <- system(paste0("bash ", tests.path, "test_cox.sh 2> /dev/null"),
              intern=TRUE)
cat("OK\n")

dose.add.PA <- read.table("coxph_dose_add.out.txt",
                          head=TRUE)[, colsAddDose]
prob.add.PA <- read.table("coxph_prob_add.out.txt",
                          head=TRUE)[, colsAddProb]
prob.dom.PA <- read.table("coxph_prob_domin.out.txt",
                          head=TRUE)[, colsDom]
prob.rec.PA <- read.table("coxph_prob_recess.out.txt",
                          head=TRUE)[, colsRec]
prob.odom.PA <- read.table("coxph_prob_over_domin.out.txt",
                           head=TRUE)[, colsOdom]
prob.2df.PA <- read.table("coxph_prob_2df.out.txt",
                          head=TRUE)[, cols2df]

####
## run analysis in R
####
attach(pheno)

cat("Comparing R output with ProbABEL output\t\t")

source(paste0(srcdir, "run_model_coxph.R"))

model.fn.0 <-
    "coxph( Surv(fupt_chd, chd)[noNA] ~ mu[noNA] + sex[noNA] + age[noNA] + othercov[noNA] )"
model.fn <- "coxph( Surv(fupt_chd, chd) ~ mu + sex + age + othercov + snp1 )"

## Additive model, dosages
snpdose <- "dose[, i]"
dose.add.R <- run.model(model.fn.0, model.fn, snpdose)
colnames(dose.add.R) <- colsAddDose
rownames(dose.add.R) <- NULL
stopifnot( all.equal(dose.add.PA, dose.add.R, tol=tol) )
cat("additive ")


## Additive model, probabilities
snpprob <- "doseFromProb[, i]"
prob.add.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.add.R) <- colsAddProb
rownames(prob.add.R) <- NULL
stopifnot( all.equal(prob.add.PA, prob.add.R, tol=tol) )
cat("additive ")

## dominant model
snpprob <- "prob[, indexHom] + prob[, indexHet]"
prob.dom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.dom.R) <- colsDom
rownames(prob.dom.R) <- NULL
stopifnot( all.equal(prob.dom.PA, prob.dom.R, tol=tol) )
cat("dominant ")

## recessive model
snpprob <- "prob[, indexHom]"
prob.rec.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.rec.R) <- colsRec
rownames(prob.rec.R) <- NULL
stopifnot( all.equal(prob.rec.PA, prob.rec.R, tol=tol) )
cat("recessive ")

## over-dominant model
snpprob <- "prob[, indexHet]"
prob.odom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.odom.R) <- colsOdom
rownames(prob.odom.R) <- NULL
stopifnot( all.equal(prob.odom.PA, prob.odom.R, tol=tol) )
cat("overdominant ")


## 2df model
model.fn <-
    "coxph( Surv(fupt_chd, chd) ~ sex + age + othercov + snp1 + snp2 )"
snpd1 <- "prob[, indexHet]"
snpd2 <- "prob[, indexHom]"
prob.2df.R <- run.model(model.fn.0, model.fn, snpd1, snpd2)
colnames(prob.2df.R) <- cols2df
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(prob.2df.PA, prob.2df.R, tol=tol) )
cat("2df\n")



cat("\nRun the same checks, but without any covariates (see bug #1266)\n")
cat("Running ProbABEL...\t\t\t\t")
tmp <- system(paste0("bash ", tests.path, "test_cox_nocovar.sh 2> /dev/null"),
              intern=TRUE)
cat("OK\n")

dose.add.PA <- read.table("coxph_dose_nocovar_add.out.txt",
                          head=TRUE)[, colsAddDose]
prob.add.PA <- read.table("coxph_prob_nocovar_add.out.txt",
                          head=TRUE)[, colsAddProb]
prob.dom.PA <- read.table("coxph_prob_nocovar_domin.out.txt",
                          head=TRUE)[, colsDom]
prob.rec.PA <- read.table("coxph_prob_nocovar_recess.out.txt",
                          head=TRUE)[, colsRec]
prob.odom.PA <- read.table("coxph_prob_nocovar_over_domin.out.txt",
                           head=TRUE)[, colsOdom]
prob.2df.PA <- read.table("coxph_prob_nocovar_2df.out.txt",
                          head=TRUE)[, cols2df]

model.fn.0 <-
    "coxph( Surv(fupt_chd, chd)[noNA] ~  mu[noNA] )"
model.fn <- "coxph( Surv(fupt_chd, chd) ~  mu + snp1 )"


cat("Comparing R output with ProbABEL output\t\t")

## Additive model, dosages
snpdose <- "dose[, i]"
dose.add.R <- run.model(model.fn.0, model.fn, snpdose)
colnames(dose.add.R) <- colsAddDose
rownames(dose.add.R) <- NULL
stopifnot( all.equal(dose.add.PA, dose.add.R, tol=tol) )
cat("additive ")


## Additive model, probabilities
snpprob <- "doseFromProb[, i]"
prob.add.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.add.R) <- colsAddProb
rownames(prob.add.R) <- NULL
stopifnot( all.equal(prob.add.PA, prob.add.R, tol=tol) )
cat("additive ")

## dominant model
snpprob <- "prob[, indexHom] + prob[, indexHet]"
prob.dom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.dom.R) <- colsDom
rownames(prob.dom.R) <- NULL
stopifnot( all.equal(prob.dom.PA, prob.dom.R, tol=tol) )
cat("dominant ")

## recessive model
snpprob <- "prob[, indexHom]"
prob.rec.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.rec.R) <- colsRec
rownames(prob.rec.R) <- NULL
stopifnot( all.equal(prob.rec.PA, prob.rec.R, tol=tol) )
cat("recessive ")

## over-dominant model
snpprob <- "prob[, indexHet]"
prob.odom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.odom.R) <- colsOdom
rownames(prob.odom.R) <- NULL
stopifnot( all.equal(prob.odom.PA, prob.odom.R, tol=tol) )
cat("overdominant ")


## 2df model
model.fn <-
    "coxph( Surv(fupt_chd, chd) ~ mu+ snp1 + snp2 )"
snpd1 <- "prob[, indexHet]"
snpd2 <- "prob[, indexHom]"
prob.2df.R <- run.model(model.fn.0, model.fn, snpd1, snpd2)
colnames(prob.2df.R) <- cols2df
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(prob.2df.PA, prob.2df.R, tol=tol) )
cat("2df\n")

cat("\t\t\t\t\t\tOK\n")
