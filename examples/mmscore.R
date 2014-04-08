#=====================================================================================
#
#       Filename:  mmscore.R
#
#    Description:  Example how to get inverse of the
#                  variance-covariance matrix from GenABEL and right
#                  phenotype table to use it in ProbABEL.
#
#        Version:  1.0
#        Created:  26_Jan-2009
#       Revision:  2010.04.18 (YA)
#
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
#          Email:  m.struchalin@@erasmusmc.nl
#
#=====================================================================================

## You have to have the GenABEL package installed on your computer
library(GenABEL)

## Load example data. Use your data here instead of example. All
## phenotypes you need must be there
data(ge03d2.clean)
data <- ge03d2.clean

## Choose snps which we are going to use as example. Just change the
## snps array if you'd like to use other snps (or use all)
snps <- c("rs70099", "rs735579", "rs9088604", "rs1413801", "rs4911638")

data <- data[!is.na(data@phdata$height),]
data <- data[!is.na(data@phdata$sex),]
data <- data[!is.na(data@phdata$age),]


## Take only 500 people for this exercise
data <- data[1:500,]

## Calculate the kinship matrix
gkin <- ibs(data[, autosomal(data)], w="freq")

## Estimate the polygenic model
h2ht <- polygenic(height~sex+age,
                  data=data,
                  kin=gkin,
                  steptol=1.e-9,
                  gradtol=1.e-9)

## Get the inverse of the variance-covariance matrix
InvSigma <- h2ht$InvSigma

## Get the phenotypes for analysis.
pheno <- data@phdata[,c("id", "height", "sex","age")]


#get rid of na
#pheno_no_na <- na.omit(pheno)

#give row names to inverse of the variance-covariance matrix
#rownames(InvSigma) <- pheno_no_na$id

## Save the inverse variance-covariance matrix it to a file. We'll use
## it in ProbABEL for mmscore
write.table(InvSigma, file="mmscore_InvSigma_aj.sex.age.dat",
            row.names=TRUE,
            col.names=FALSE,
            quote=FALSE)

## Get residuals from analysis, based on covariate effects only.
height_residuals <- h2ht$residualY

## Create a table with two columns: id and trait
pheno_residuals <- data.frame(id=pheno$id,
                              height_residuals=height_residuals)

## Add row names
rownames(pheno_residuals) <- as.character(pheno_residuals$id)

## Save it into the file. We will use this file in ProbABEL
write.table(pheno_residuals,
            file="mmscore_pheno.PHE",
            row.names=FALSE,
            quote=FALSE)

## Now we have two files:
## 1) inverse of the variance-covariance matrix
## 2) residuals of the phenotype, which will be the new phenotype that
## ProbABEL will analyse.

## Now, go to ProbABEL and start analysis

