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
tmp <- system(paste0("cd ", tests.path, "; bash test_bt.sh; cd -"),
              intern=TRUE)
cat("OK\n")

resPaAddDose <- read.table(
    paste0(tests.path, "logist_add.out.txt"),
    head=TRUE)[,
        c("beta_SNP_add",
          "sebeta_SNP_add",
          "chi2_SNP")]
resPaAddProb <- read.table(
    paste0(tests.path, "logist_prob_add.out.txt"),
    head=TRUE)[, c("beta_SNP_addA1",
        "sebeta_SNP_addA1",
        "chi2_SNP_A1")]
resPaDom <- read.table(
    paste0(tests.path, "logist_prob_domin.out.txt"),
    head=TRUE)[,
        c("beta_SNP_domA1",
          "sebeta_SNP_domA1",
          "chi2_SNP_domA1")]
resPaRec <- read.table(
    paste0(tests.path, "logist_prob_recess.out.txt"),
    head=TRUE)[,
        c("beta_SNP_recA1",
          "sebeta_SNP_recA1",
          "chi2_SNP_recA1")]
resPaOdom <- read.table(
    paste0(tests.path, "logist_prob_over_domin.out.txt"),
    head=TRUE)[,
        c("beta_SNP_odomA1",
          "sebeta_SNP_odomA1",
          "chi2_SNP_odomA1")]
resPa2df <- read.table(
    paste0(tests.path, "logist_prob_2df.out.txt"),
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
## run analysis on dose
dose.add.R <- data.frame(beta_SNP_add   = numeric(),
                         sebeta_SNP_add = numeric(),
                         chi2_SNP       = numeric())
prob.add.R <- data.frame(beta_SNP_addA1   = numeric(),
                         sebeta_SNP_addA1 = numeric(),
                         chi2_SNP_A1      = numeric())
for (i in 3:dim(dose)[2]) {
        noNA <- which( !is.na(dose[, i]) )
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial)

        model.dose <- glm( chd ~ sex + age + othercov + dose[, i],
                          family=binomial)
        sm.dose <- summary(model.dose)$coef[5, 1:2]

        model.prob <- glm( chd ~ sex + age + othercov + doseFromProb[,
                                                                     i],
                          family=binomial )
        sm.prob <- summary(model.prob)$coef[5, 1:2]

        lrt.dose <- 2 * ( logLik( model.dose ) - logLik( model.0 ) )
        lrt.prob <- 2 * ( logLik( model.prob ) - logLik( model.0 ) )

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
        regProb  <- prob[, indexHom] + prob[, indexHet]

        noNA    <- which( !is.na(regProb) )
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial )
        model   <- glm( chd ~ sex + age + othercov + regProb,
                       family=binomial )
        sm      <- summary(model)$coef[5, 1:2]
        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

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
        regProb <- prob[, indexHom]

        noNA    <- which( !is.na(regProb) )
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial )
        model   <- glm( chd ~ sex + age + othercov + regProb,
                       family=binomial )
        sm      <- summary(model)$coef[5, 1:2]
        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

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
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial )
        model   <- glm( chd ~ sex + age + othercov + regProb,
                       family=binomial )
        sm      <- summary(model)$coef[5, 1:2]
        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

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
        model.0 <- glm( chd[noNA] ~ sex[noNA] + age[noNA] + othercov[noNA],
                       family=binomial )
        model   <- glm( chd ~ sex + age + othercov + prob[, indexHet] +
                       prob[, indexHom], family=binomial )
        smA1A2  <- summary(model)$coef[5, 1:2]
        smA1A1  <- summary(model)$coef[6, 1:2]
        lrt     <- 2 * ( logLik( model ) - logLik( model.0 ) )

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
