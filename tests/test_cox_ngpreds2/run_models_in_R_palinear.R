###
# load the data
###

# load phenotypic data
pheno <- read.table("../../examples/height.txt", head=TRUE, string=FALSE)

# load genetic DOSE data
dose <- read.table("../../examples/test.mldose.2", head=FALSE, string=FALSE)
# remove "1->" from the names of dose-IDs
idNames <- dose[, 1]
idNames <- sub("[0-9]+->", "", idNames)
dose[, 1] <- idNames
# check consistency of names
table(dose[, 1] == pheno[, 1])

# load genetic PROB data
prob <- read.table("../../examples/test.mlprob", head=FALSE, string=FALSE)
# remove "1->" from the names of prob-IDs
idNames <- prob[, 1]
idNames <- sub("[0-9]+->", "", idNames)
prob[, 1] <- idNames
# check consistency of names
table(prob[, 1] == pheno[, 1])

# check consistency DOSE <-> PROB
doseFromProb <- matrix(NA, ncol=dim(dose)[2], nrow=dim(dose)[1])
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        doseFromProb[, i] <- prob[, indexHom] * 2 + prob[, indexHet]
        print( table( ( dose[, i] - doseFromProb[, i] ) > 1e-8 ) )
}

###
# run analysis
###

attach(pheno)

# run analysis on dose
for (i in 3:dim(dose)[2]) {
        smr <- summary( lm( height ~ sex + age + dose[, i] ) )
        print( smr$coef[4, 1:2] )
        smr <- summary( lm( height ~ sex + age + doseFromProb[, i] ) )
        print( smr$coef[4, 1:2] )
}

# dominant model
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHom] + prob[, indexHet]
        smr <- summary( lm( height ~ sex + age + regProb ) )
        print( smr$coef[4, 1:2] )
}

# recessive model
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHom]
        smr <- summary( lm( height ~ sex + age + regProb ) )
        print( smr$coef[4, 1:2] )
}

# over-dominant model
for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        regProb <- prob[, indexHet]
        smr <- summary( lm( height ~ sex + age + regProb ) )
        print( smr$coef[4, 1:2] )
}
