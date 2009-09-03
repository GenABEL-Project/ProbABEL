#=====================================================================================
#
#       Filename:  mmscore.R
#
#    Description:  Example how to get inverse of the variance-covariance matrix from GenABEL and right phenotype table to use it in ProbABEL.
#
#        Version:  1.0
#        Created:  26_Jan-2009
#       Revision:  none
#       
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
#          Email:  m.struchalin@@erasmusmc.nl
#
#=====================================================================================


#You have to have GenABEL installed on your computer
library(GenABEL)
data(ge03d2.clean)





#load example data. Use your data here instead of example. All phenotypes you need must be there
data(ge03d2.clean)


data <- ge03d2.clean



#choose snps which we are going to use as example. Just change snps array if you'd like to use other snps
snps <- c("rs70099","rs735579","rs9088604","rs1413801","rs4911638")


data <- data[!is.na(data@phdata$height),]
data <- data[!is.na(data@phdata$sex),]
data <- data[!is.na(data@phdata$age),]



for (snp in snps) {
   data <- data[!is.na(as.numeric(data@gtdata[,snp])),]
   }

data <- data[1:500,]



#calculate kinship matrix
gkin <- ibs(data[,autosomal(data)],w="freq")









#estimate polygenic model
h2ht <- polygenic(height~sex+age,data=data,kin=gkin, steptol=1.e-9,gradtol=1.e-9)


#get inverse of the variance-covariance matrix
InvSigma <- h2ht$InvSigma


#get phenotypes for analysis.
pheno <- data@phdata[,c("id", "height","sex","age")]


#get rid of na
pheno_no_na <- na.omit(pheno)

#give row names to inverse of the variance-covariance matrix
rownames(InvSigma) <- pheno_no_na$id


#save it to file. We'll use it in ProbABEL for mmscore
write.table(InvSigma, file="mmscore_InvSigma_aj.sex.age.dat", col.names=F, quote=F)


#Get residuals from analysis, based on covariate effects only.
height_residuals <- h2ht$residualY

#Create table with two columns: id and trait  
pheno_residuals <- data.frame(id=pheno$id, height_residuals=height_residuals)


#give row names 
rownames(pheno_residuals) <- as.character(pheno_residuals$id)



#save it into the file. We will use this file in ProbABEL
write.table(pheno_residuals, file="mmscore_pheno.PHE", row.names=F, quote=F)




#Now we have two files 
#1) inverse of the variance-covariance matrix
#2) residuals

#go to ProbABEL and start analysis






#____________________________________________________

#Create test file with genotypes

data_cut <- data[,snps]

gen_table <- as.numeric(data_cut)


#Replace NA by mean for each snp. NA is forbiden in genotypes in ProbABEL input (!).
for(snpnum in 1:dim(gen_table)[2]) 
	{
	mean <- mean(gen_table[,snpnum], na.rm=T)	
	gen_table[is.na(gen_table[,snpnum]),snpnum]	<- mean
	}


gen_table_df <- data.frame(gen_table)


gen_table_df$MLDOSE <- gen_table_df[,1]
gen_table_df[,1] <- "MLDOSE"

colnam <- colnames(gen_table_df)


#colnam[-c(1:length(colnam)-1)] <- colnam[1]
#colnam[1] <- "MLDOSE"

colnam[length(colnam)] <- colnam[1]
colnam[1] <- "MLDOSE"


colnames(gen_table_df) <- colnam 


rownames <- rownames(gen_table_df)
rownames <- paste("1->", rownames, sep="")
rownames(gen_table_df) <- rownames

write.table(file="mmscore_gen.mldose", gen_table_df, col.names=F, quote=F)

mlinfo <- data.frame(SNP=colnam[2:length(colnam)])
mlinfo$Al1 <- "A"
mlinfo$Al2 <- "B"
mlinfo$Freq1 <- 0.5847
mlinfo$MAF <- 0.5847
mlinfo$Quality <- 0.5847
mlinfo$Rsq <- 0.5847

write.table(mlinfo, "mmscore_gen.mlinfo", row.names=F, quote=F)

