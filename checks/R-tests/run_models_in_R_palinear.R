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
prnt("Running ProbABEL...")
tmp <- system(paste0("bash ", tests.path, "test_qt.sh"),
              intern=TRUE)
cat("OK\n")

dose.add.PA <- read.table("linear_base_add.out.txt",
                          head=TRUE)[, colsAdd]
prob.add.PA <- read.table("linear_ngp2_add.out.txt",
                          head=TRUE)[, colsAdd]
prob.dom.PA <- read.table("linear_ngp2_domin.out.txt",
                          head=TRUE)[, colsDom]
prob.rec.PA <- read.table("linear_ngp2_recess.out.txt",
                          head=TRUE)[, colsRec]
prob.odom.PA <- read.table("linear_ngp2_over_domin.out.txt",
                           head=TRUE)[, colsOdom]
prob.2df.PA <- read.table("linear_ngp2_2df.out.txt",
                          head=TRUE)[, cols2df]

## Fix betas, sebetas, chi^2 for the case that there is no variation
## (SNP 6 in the info file). ProbABEL lists them all as 0.0, R lists
## them as:
prob.dom.PA[6, 2:4] <- c(NaN, NaN, 0.0)

## For the 2df model the last SNP is interchangeable: ProbABEL
## calculates the beta for the other SNP than R. This causes the beta
## to have the wrong sign. This part changes the position of the SNP
## beta (and swaps sign) and SE if beta and the other SE are 0.
if (sum(abs(prob.2df.PA[6, 2:3])) == 0){
    prob.2df.PA[6, 2:3] <- c(prob.2df.PA[6, 4] * -1, prob.2df.PA[6, 5])
    prob.2df.PA[6, 4:5] <- c(NA, NA)
}

####
## run analysis in R
####
attach(pheno)

cat("Comparing R output with ProbABEL output:\n")

source(paste0(srcdir, "run_model_linear.R"))

model.fn.0 <- "lm( height[noNA] ~ sex[noNA] + age[noNA] )"
model.fn   <- "lm( height ~ sex + age + snp )"

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
#print("OK", justify="right")
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
prob.2df.R <- data.frame()
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- lm( height[noNA] ~ sex[noNA] + age[noNA] )
        model   <- lm( height ~ sex + age + prob[, indexHet] + prob[, indexHom] )
        smA1A2  <- summary(model)$coef[4, 1:2]

        ## When all coefficients are NA, they don't show up in $coeff
        if ( nrow(summary(model)$coeff) > 4 ) {
            smA1A1 <- summary(model)$coef[5, 1:2]
        } else {
            smA1A1 <- c(NA, NA)
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
cat("OK\n")
