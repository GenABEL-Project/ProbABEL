#
# file containing all phenotypes
# variables should be separated by space or tab, missing values coded as NA
# header line should contain variable names. 
# Compulsory key variable should be named "id"
# All people listed in ids_order_file should be present in original_phenofile
#
original_phenofile <- "allheight.txt"

#
# IDs present in mldose file (in mldoes file order)
# normally generated with getIDS.pl
#
ids_order_file <- "mldose.IDS"

#
# output file names (the one used by linear/logistic)
#
output_phenofile <- "height.txt"

#
# trait to be analysied 
# (name should be the same as in the header line of original_phenofile)
# the analysed trait comes first
#
trait <- c("height")

#
# covariates to be included into model
# (names should be the same as in the header line of original_phenofile)
#
covars <- c("sex","age")

#
# code
#

ophe <- read.table(original_phenofile,head=T,strings=F)
rownames(ophe) <- ophe$id
idso <- scan(ids_order_file,what=character())

newphe <- ophe[idso,c("id",trait,covars)]

write.table(newphe,file=output_phenofile)

### convert mach 2 fvf

library(GenABEL)
mach2databel("test.mldose","test.mlinfo","test.mldose_fvf")
mach2databel("mmscore_gen.mldose","mmscore_gen.mlinfo","mmscore_gen.mldose_fvf")

