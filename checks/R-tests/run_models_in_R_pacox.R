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
prnt("Running ProbABEL...")
tmp <- system(paste0("bash ", tests.path, "test_cox.sh 2> /dev/null"),
              intern=TRUE)
cat("OK\n")

dose.add.PA <- read.table("coxph_dose_add.out.txt",
                          head=TRUE)[, colsAdd]
prob.add.PA <- read.table("coxph_prob_add.out.txt",
                          head=TRUE)[, colsAdd]
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

cat("Comparing R output with ProbABEL output:\n")

source(paste0(srcdir, "run_model_coxph.R"))

model.fn.0 <-
    "coxph( Surv(fupt_chd, chd)[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA] )"
model.fn <- "coxph( Surv(fupt_chd, chd) ~ sex + age + othercov + snp1 )"

## Additive model, dosages
prnt(" additive (dosages)")
snpdose <- "dose[, i]"
dose.add.R <- run.model(model.fn.0, model.fn, snpdose)
colnames(dose.add.R) <- colsAdd
rownames(dose.add.R) <- NULL
stopifnot( all.equal(dose.add.PA, dose.add.R, tol=tol) )
cat("OK\n")

## Additive model, probabilities
prnt(" additive (probabilities)")
snpprob <- "doseFromProb[, i]"
prob.add.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.add.R) <- colsAdd
rownames(prob.add.R) <- NULL
stopifnot( all.equal(prob.add.PA, prob.add.R, tol=tol) )
cat("OK\n")

## dominant model
prnt(" dominant")
snpprob <- "prob[, indexHom] + prob[, indexHet]"
prob.dom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.dom.R) <- colsDom
rownames(prob.dom.R) <- NULL
stopifnot( all.equal(prob.dom.PA, prob.dom.R, tol=tol) )
cat("OK\n")

## recessive model
prnt(" recessive")
snpprob <- "prob[, indexHom]"
prob.rec.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.rec.R) <- colsRec
rownames(prob.rec.R) <- NULL
stopifnot( all.equal(prob.rec.PA, prob.rec.R, tol=tol) )
cat("OK\n")

## over-dominant model
prnt(" overdominant")
snpprob <- "prob[, indexHet]"
prob.odom.R <- run.model(model.fn.0, model.fn, snpprob)
colnames(prob.odom.R) <- colsOdom
rownames(prob.odom.R) <- NULL
stopifnot( all.equal(prob.odom.PA, prob.odom.R, tol=tol) )
cat("OK\n")

## 2df model
prnt(" 2df")
model.fn <-
    "coxph( Surv(fupt_chd, chd) ~ sex + age + othercov + snp1 + snp2 )"
snpd1 <- "prob[, indexHet]"
snpd2 <- "prob[, indexHom]"
prob.2df.R <- run.model(model.fn.0, model.fn, snpd1, snpd2)
colnames(prob.2df.R) <- cols2df
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(prob.2df.PA, prob.2df.R, tol=tol) )
cat("OK\n")
