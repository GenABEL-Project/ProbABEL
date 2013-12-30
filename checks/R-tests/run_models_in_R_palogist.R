cat("Checking logistic regression...\n")

args <- commandArgs(TRUE)
srcdir <- args[1]

if (is.na(srcdir)) {
    srcdir <- "./"
}

pheno.file <- "logist_data.txt"

source(paste0(srcdir, "initial_checks.R"))

####
## Run ProbABEL to get the output data we want to compare/verify
####
cat("Running ProbABEL...\t\t\t\t")
tmp <- system(paste0("bash ", tests.path, "test_bt.sh"),
              intern=TRUE)
cat("OK\n")

dose.add.PA <- read.table("logist_add.out.txt",
                          head=TRUE)[, colsAddDose]
prob.add.PA <- read.table("logist_prob_add.out.txt",
                          head=TRUE)[, colsAddProb]
prob.dom.PA <- read.table("logist_prob_domin.out.txt",
                          head=TRUE)[, colsDom]
prob.rec.PA <- read.table("logist_prob_recess.out.txt",
                          head=TRUE)[, colsRec]
prob.odom.PA <- read.table("logist_prob_over_domin.out.txt",
                           head=TRUE)[, colsOdom]
prob.2df.PA <- read.table("logist_prob_2df.out.txt",
                          head=TRUE)[, cols2df]

## Fix chi^2 for the case that there is no variation (SNP 6 in the
## info file). ProbABEL lists it as NaN, R lists it as:
prob.dom.PA[6, 4] <- 0.0

####
## run analysis in R
####
attach(pheno)

cat("Comparing R output with ProbABEL output\t\t")

source(paste0(srcdir, "run_model_logist.R"))

model.fn.0 <-
    "glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA], family=binomial)"
model.fn  <- "glm( chd ~ sex + age + othercov + snp, family=binomial )"

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
prob.2df.R <- data.frame()
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial )
        model   <- glm( chd ~ sex + age + othercov + prob[, indexHet] +
                       prob[, indexHom], family=binomial )
        smA1A2  <- summary(model)$coef[5, 1:2]

        if ( nrow(summary(model)$coeff) > 5 ) {
            smA1A1 <- summary(model)$coef[6, 1:2]
        } else {
            smA1A1 <- c(NaN, NaN)
        }

        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

        rsq <- prob.2df.PA[i-2, "Rsq"]
        if( rsq < rsq.thresh) {
            row <- c(rsq, NaN, NaN, NaN, NaN, NaN)
        } else {
            row <- c(rsq, smA1A2[1], smA1A2[2], smA1A1[1], smA1A1[2], lrt)

        }
        prob.2df.R <- rbind(prob.2df.R, row)
}
colnames(prob.2df.R) <- cols2df
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(prob.2df.PA, prob.2df.R, tol=tol) )
cat("2df\n")

cat("\t\t\t\t\t\tOK\n")
