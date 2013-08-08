library(survival)

###
# load the data
###

# load phenotypic data
pheno <- read.table("../../examples/coxph_data.txt", head=TRUE,
                    string=FALSE)

# load genetic DOSE data
dose <- read.table("../../examples/test.mldose", head=FALSE,
                   string=FALSE)
# remove "1->" from the names of dose-IDs
idNames <- dose[, 1]
idNames <- sub("[0-9]+->", "", idNames)
dose[,1] <- idNames
# check consistency of names
table(dose[, 1] == pheno[, 1])

# load genetic PROB data
prob <- read.table("../../examples/test.mlprob", head=FALSE, string=FALSE)
# remove "1->" from the names of prob-IDs
idNames <- prob[, 1]
idNames <- sub("[0-9]+->", "", idNames)
prob[,1] <- idNames
# check consistency of names
table(prob[,1] == pheno[,1])

# check consistency DOSE <-> PROB
doseFromProb <- matrix(NA, ncol=dim(dose)[2], nrow=dim(dose)[1])
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        doseFromProb[, i] <- prob[, indexHom] * 2 + prob[, indexHet]
        print( table( ( dose[,i] - doseFromProb[,i] ) > 1e-8 ) )
}

###
# run analysis
###

# run ProbABEL
system("cd ../../examples/; sh example_cox.sh; cd -")
resPaAddDose <-
  read.table("../../examples/coxph_dose_add.out.txt",
             head=TRUE)[, c("beta_SNP_add","sebeta_SNP_add")]
resPaAddProb <-
  read.table("../../examples/coxph_prob_add.out.txt",
             head=TRUE)[, c("beta_SNP_addA1","sebeta_SNP_addA1")]
resPaDom <-
  read.table("../../examples/coxph_prob_domin.out.txt",
             head=TRUE)[, c("beta_SNP_domA1","sebeta_SNP_domA1")]
resPaRec <-
  read.table("../../examples/coxph_prob_recess.out.txt",
             head=TRUE)[, c("beta_SNP_recA1","sebeta_SNP_recA1")]
resPaOdom <-
  read.table("../../examples/coxph_prob_over_domin.out.txt",
             head=TRUE)[, c("beta_SNP_odom","sebeta_SNP_odom")]

attach(pheno)

# run analysis on dose
emptyRes <- resPaAddDose
emptyRes[] <- NA
resRAddDose <- resRAddProb <- resRDom <- resRRec <- resROdom <- emptyRes
for (i in 3:dim(dose)[2]) {
        smr <- summary( coxph(Surv(fupt_chd, chd) ~ sex + age +
                              othercov + dose[,i] ) )
        resRAddDose[i-2,] <-  smr$coef[4, c(1,3)]
        smr <- summary( coxph(Surv(fupt_chd, chd) ~ sex + age +
                              othercov + doseFromProb[,i] ) )
        resRAddProb[i-2,] <- smr$coef[4, c(1,3)]
}
cat("additive model (dose):\n")
resPaAddDose/resRAddDose
cat("additive model (prob):\n")
resPaAddProb/resRAddProb

cat("dominant model (prob):\n")
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[,indexHom] + prob[,indexHet]
        smr <- summary( coxph(Surv(fupt_chd, chd) ~ sex + age +
                              othercov + regProb ) )
        resRDom[i-2,] <- smr$coef[4, c(1,3)]
}
resPaDom/resRDom

cat("recessive model (prob):\n")
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[,indexHom]
        smr <- summary( coxph(Surv(fupt_chd, chd) ~ sex + age +
                              othercov + regProb ) )
        resRRec[i-2,] <- smr$coef[4, c(1,3)]
}
resPaRec/resRRec

cat("over-dominant model (prob):\n")
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[,indexHet]
        smr <- summary( coxph(Surv(fupt_chd, chd) ~ sex + age +
                              othercov + regProb ) )
        resROdom[i-2,] <- smr$coef[4, c(1,3)]
}
resPaOdom/resROdom
