# This R script serves as an example of how to create input files for
# running ProbABEL.
# It requires several files (listed below) and basically selects the
# phenotype and covariates you want from a larger phenotype file,
# orders them in the same order as the mldose file. At the end it
# converts existing dose and prob files to DatABEL format.

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
# normally generated with extIDS.pl
#
ids_order_file <- "mldose.IDS"

#
# output file names (the one used by linear/logistic)
#
output_phenofile <- "height.txt"

#
# trait to be analysed
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

ophe <- read.table(original_phenofile, header=TRUE, stringsAsFactors=FALSE)
rownames(ophe) <- ophe$id
idso <- scan(ids_order_file, what=character())

newphe <- ophe[idso, c("id", trait,covars)]

write.table(newphe, file=output_phenofile, quote=FALSE, row.names=FALSE)

if (!require(GenABEL))
  stop("further code requires the GenABEL library to be installed")
if (!require(DatABEL))
  stop("further code requires the DatABEL library to be installed")


fvdose <- mach2databel(imputedg="test.mldose",
                       mlinfo="test.mlinfo",
                       out="test.dose")
fvprob <- mach2databel(imputedg="test.mlprob",
                       mlinfo="test.mlinfo",
                       out="test.prob",
                       isprob=TRUE)
mmdose <- mach2databel("mmscore_gen.mldose",
                       "mmscore_gen.mlinfo",
                       "mmscore_gen.dose")
mmprob <- mach2databel("mmscore_gen.mlprob",
                       "mmscore_gen.mlinfo",
                       "mmscore_gen.prob",
                       isprob=TRUE)
