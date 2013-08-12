cat("Checking Cox PH regression...\n")
library(survival)

## Set tolerance for comparing various outputs
tol <- 1e-5

###
# load the data
###
example.path <- "../../examples/"

# load phenotypic data
pheno <- read.table("../../examples/coxph_data.txt", head=TRUE,
                    string=FALSE)

# load genetic DOSE data
dose <- read.table("../../examples/test.mldose",
                   head=FALSE, string=FALSE)
# remove "1->" from the names of dose-IDs
idNames   <- dose[, 1]
idNames   <- sub("[0-9]+->", "", idNames)
dose[, 1] <- idNames
cat("Dose: check consistency of names\t\t")
stopifnot( all.equal(dose[, 1], pheno[, 1], tol) )
cat("OK\n")

# load genetic PROB data
prob <- read.table("../../examples/test.mlprob",
                   head=FALSE, string=FALSE)
# remove "1->" from the names of prob-IDs
idNames   <- prob[, 1]
idNames   <- sub("[0-9]+->", "", idNames)
prob[, 1] <- idNames
cat("Prob: check consistency of names\t\t")
stopifnot( all.equal(prob[, 1], pheno[, 1], tol) )
cat("OK\n")

# check consistency DOSE <-> PROB
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

####
## Run ProbABEL to get the output data we want to compare/verify
####
cat("Running ProbABEL...\t\t\t\t")
tmp <- system("cd ../../examples/; sh example_cox.sh; cd -",
              intern=TRUE)
cat("OK\n")

resPaAddDose <- read.table(
    paste0(example.path, "coxph_dose_add.out.txt"),
    head=TRUE)[,
        c("beta_SNP_add",
          "sebeta_SNP_add",
          "chi2_SNP")]
resPaAddProb <- read.table(
    paste0(example.path, "coxph_prob_add.out.txt"),
    head=TRUE)[, c("beta_SNP_addA1",
        "sebeta_SNP_addA1",
        "chi2_SNP_A1")]
resPaDom <- read.table(
    paste0(example.path, "coxph_prob_domin.out.txt"),
    head=TRUE)[,
        c("beta_SNP_domA1",
          "sebeta_SNP_domA1",
          "chi2_SNP_domA1")]
resPaRec <- read.table(
    paste0(example.path, "coxph_prob_recess.out.txt"),
    head=TRUE)[,
        c("beta_SNP_recA1",
          "sebeta_SNP_recA1",
          "chi2_SNP_recA1")]
resPaOdom <- read.table(
    paste0(example.path, "coxph_prob_over_domin.out.txt"),
    head=TRUE)[,
        c("beta_SNP_odomA1",
          "sebeta_SNP_odomA1",
          "chi2_SNP_odomA1")]
resPa2df <- read.table(
    paste0(example.path, "coxph_prob_2df.out.txt"),
    head=TRUE)[,
        c("beta_SNP_A1A2",
          "sebeta_SNP_A1A2",
          "beta_SNP_A1A1",
          "sebeta_SNP_A1A1",
          "chi2_SNP_2df")]

####
## run analysis in R
####
attach(pheno)

cat("Comparing R output with ProbABEL output\t\t")
# run analysis on dose
## emptyRes <- resPaAddDose
## emptyRes[] <- NA
## resRAddDose <- resRAddProb <- resRDom <- resRRec <- resROdom <- emptyRes

dose.add.R <- data.frame(beta_SNP_add   = numeric(),
                         sebeta_SNP_add = numeric(),
                         chi2_SNP       = numeric())
prob.add.R <- data.frame(beta_SNP_addA1   = numeric(),
                         sebeta_SNP_addA1 = numeric(),
                         chi2_SNP_A1      = numeric())
for (i in 3:dim(dose)[2]) {
        noNA <- which( !is.na(dose[, i]) )
        model.0 <- coxph(
            Surv(fupt_chd, chd)[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA] )

        model.dose <- coxph(Surv(fupt_chd, chd) ~ sex + age +
                            othercov + dose[,i] )
        sm.dose <- summary(model.dose)$coef[4, c(1,3)]

        model.prob <- coxph(Surv(fupt_chd, chd) ~ sex + age +
                            othercov + doseFromProb[,i] )
        sm.prob <- summary(model.prob)$coef[4, c(1,3)]

        lrt.dose <- 2 * ( model.dose$loglik[2] - model.0$loglik[2] )
        lrt.prob <- 2 * ( model.prob$loglik[2] - model.0$loglik[2] )

        row <- data.frame(
            beta_SNP_add = sm.dose[1],
            sebeta_SNP_add = sm.dose[2],
            chi2_SNP = lrt.dose)
        dose.add.R <- rbind(dose.add.R, row)

        row <- data.frame(
            beta_SNP_addA1 = sm.prob[1],
            sebeta_SNP_addA1 = sm.prob[2],
            chi2_SNP_A1 = lrt.prob)
        prob.add.R <- rbind(prob.add.R, row)
}
rownames(dose.add.R) <- NULL
rownames(prob.add.R) <- NULL
stopifnot( all.equal(resPaAddDose, dose.add.R, tol=tol) )

stopifnot( all.equal(resPaAddProb, prob.add.R, tol=tol) )


## dominant model
prob.dom.R <- data.frame(beta_SNP_domA1   = numeric(),
                         sebeta_SNP_domA1 = numeric(),
                         chi2_SNP_domA1   = numeric())
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[,indexHom] + prob[,indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- coxph( Surv(fupt_chd, chd)[noNA] ~ sex[noNA] +
                         age[noNA] + othercov[noNA])
        model   <- coxph( Surv(fupt_chd, chd) ~ sex + age +
                         othercov + regProb )
        sm      <- summary(model)$coef[4, c(1,3)]
        lrt     <- 2 * ( model$loglik[2] - model.0$loglik[2] )

        row <- data.frame(
            beta_SNP_domA1   = sm[1],
            sebeta_SNP_domA1 = sm[2],
            chi2_SNP_domA1   = lrt)
        prob.dom.R <- rbind(prob.dom.R, row)
}
rownames(prob.dom.R) <- NULL
stopifnot( all.equal(resPaDom, prob.dom.R, tol=tol) )

## recessive model
prob.rec.R <- data.frame(beta_SNP_recA1   = numeric(),
                         sebeta_SNP_recA1 = numeric(),
                         chi2_SNP_recA1   = numeric())
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[,indexHom]

        noNA    <- which( !is.na(regProb) )
        model.0 <- coxph( Surv(fupt_chd, chd)[noNA] ~ sex[noNA] +
                         age[noNA] + othercov[noNA])
        model   <- coxph( Surv(fupt_chd, chd) ~ sex + age +
                         othercov + regProb )
        sm      <- summary(model)$coef[4, c(1,3)]
        lrt     <- 2 * (  model$loglik[2] - model.0$loglik[2] )

        row <- data.frame(
            beta_SNP_recA1   = sm[1],
            sebeta_SNP_recA1 = sm[2],
            chi2_SNP_recA1   = lrt)
        prob.rec.R <- rbind(prob.rec.R, row)
}
rownames(prob.rec.R) <- NULL
stopifnot( all.equal(resPaRec, prob.rec.R, tol=tol) )


## over-dominant model
prob.odom.R <- data.frame(beta_SNP_odomA1   = numeric(),
                         sebeta_SNP_odomA1 = numeric(),
                         chi2_SNP_odomA1   = numeric())
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- coxph( Surv(fupt_chd, chd)[noNA] ~ sex[noNA] +
                         age[noNA] + othercov[noNA])
        model   <- coxph( Surv(fupt_chd, chd) ~ sex + age +
                         othercov + regProb )
        sm      <- summary(model)$coef[4, c(1,3)]
        lrt     <- 2 * (  model$loglik[2] - model.0$loglik[2] )

        row <- data.frame(
            beta_SNP_odomA1   = sm[1],
            sebeta_SNP_odomA1 = sm[2],
            chi2_SNP_odomA1   = lrt)
        prob.odom.R <- rbind(prob.odom.R, row)
}
rownames(prob.odom.R) <- NULL
stopifnot( all.equal(resPaOdom, prob.odom.R, tol=tol) )

## 2df model
prob.2df.R <- data.frame(beta_SNP_A1A2   = numeric(),
                          sebeta_SNP_A1A2 = numeric(),
                          beta_SNP_A1A1   = numeric(),
                          sebeta_SNP_A1A1 = numeric(),
                         chi2_SNP_2df     = numeric())
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- coxph( Surv(fupt_chd, chd)[noNA] ~ sex[noNA] +
                         age[noNA] + othercov[noNA])
        model   <- coxph( Surv(fupt_chd, chd) ~ sex + age +
                         othercov + prob[, indexHet] + prob[, indexHom] )
        smA1A2  <- summary(model)$coef[4, c(1,3)]
        smA1A1  <- summary(model)$coef[5, c(1,3)]
        lrt     <- 2 * (  model$loglik[2] - model.0$loglik[2] )

        row <- data.frame(
            beta_SNP_A1A2   = smA1A2[1],
            sebeta_SNP_A1A2 = smA1A2[2],
            beta_SNP_A1A1   = smA1A1[1],
            sebeta_SNP_A1A1 = smA1A1[2],
            chi2_SNP_2df   = lrt)
        prob.2df.R <- rbind(prob.2df.R, row)
}
rownames(prob.2df.R) <- NULL
stopifnot( all.equal(resPa2df, prob.2df.R, tol=tol) )

cat("OK\n")
