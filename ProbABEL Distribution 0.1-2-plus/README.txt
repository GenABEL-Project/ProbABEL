ProbABEL Version 0.1-2-plus
Modified by Han Chen (hanchen@bu.edu) on Nov 9, 2009 to extract the covariance between the estimate of beta(SNP) and the estimate of beta(interaction), based on ProbABEL version 0.1-2 (released on Oct 19, 2009) downloaded from http://mga.bionet.nsc.ru/~yurii/ABEL/

Attention:
For model Y = b_0 + b_cov1 * cov1 + b_cov2 * cov2 + ... + b_covN * covN + b_SNP * SNP + b_covX_SNP * covX * SNP
(1<=X<=N, covX can be any covariate in the phenotype file, from cov1 to covN, option --interaction=X, see ProbABEL manual for details)
Or model Y = b_0 + b_cov1 * cov1 + b_cov2 * cov2 + ... + b_covX-1 * covX-1 + b_covX+1 * covX+1 + ... + b_covN * covN + b_SNP * SNP + b_covX_SNP * covX * SNP
(1<=X<=N, covX can be any covariate in the phenotype file, from cov1 to covN, option --interaction_only=X, see ProbABEL manual for details)
This "plus" version reports naive covariance (default) or robust covariance (option --robust) between b_SNP and b_covX_SNP estimates, for palinear and palogist (can NOT report covariance for pacoxph, Cox regression).
Also, in order to get the covariance, users should NOT use the following options: --score, --allcov, --mmscore
LINEAR REGRESSION OR LOGISTIC REGRESSION ONLY!!!!!





Steps:
1. Extract the original source code, starting with ProbABEL_0.1-2.tar.gz (released on Oct 19, 2009) downloaded from http://mga.bionet.nsc.ru/~yurii/ABEL/
[user@server ~]$ tar -xzvf ProbABEL_0.1-2.tar.gz

2. Copy the patch file ProbABEL_0.1-2-plus.patch to the same directory.
[user@server ~]$ ls
ProbABEL_0.1-2  ProbABEL_0.1-2-plus.patch  ProbABEL_0.1-2.tar.gz

3. Change to ProbABEL_0.1-2 directory.
[user@server ~]$ cd ProbABEL_0.1-2
[user@server ~/ProbABEL_0.1-2]$

4. Patch the patch file.
[user@server ~/ProbABEL_0.1-2]$ patch -b -p1 < ../ProbABEL_0.1-2-plus.patch
patching file src/main.cpp
patching file src/reg1.h
patching file src/version.h

5. Compile.
[user@server ~/ProbABEL_0.1-2]$ make
g++ -I src/include -O3 -DLOGISTIC src/main.cpp -o bin/palogist
g++ -I src/include -O3 -DLINEAR src/main.cpp -o bin/palinear
g++ -I src/include -O3 -DCOXPH src/chinv2.c src/cholesky2.c src/chsolve2.c src/dmatrix.c src/coxfit2.c src/main.cpp -o bin/pacoxph
cp src/extIDS.pl bin/.
cp src/prepare_data.R bin/.
cp src/probabel.pl bin/probabel.pl_example
cp src/probabel_config.cfg bin/probabel_config.cfg_example





Note:
None of the calculations implemented in the original version of ProbABEL (0.1-2) is disabled in this 0.1-2-plus version.
Users should expect an additional column (covariance) in the output of all models except 2df model, if they run linear regression or logistic regression.
Two additional columns of covariances will appear in the output of 2df model, if an MLPROB file is used as the genomic predictor file (option --ngpreds=2).
Users can still use Cox regression, or score test in palinear and palogist, the output will be the same as the output of the original version (0.1-2), no covariance will be reported.

Examples:
[user@server ~/data_ex]$ ls
ex.mldose  ex.mlinfo  ex.mlprob  pheno.txt
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mldose -i ex.mlinfo --out=ex1 --interaction=1
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mlprob -i ex.mlinfo --ngpreds=2 --out=ex2 --interaction=1
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mldose -i ex.mlinfo --out=ex3 --interaction_only=1
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mlprob -i ex.mlinfo --ngpreds=2 --out=ex4 --interaction_only=1
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mldose -i ex.mlinfo --out=ex5 --interaction=1 --robust
[user@server ~/data_ex]$ ../ProbABEL_0.1-2/bin/palinear -p pheno.txt -d ex.mlprob -i ex.mlinfo --ngpreds=2 --out=ex6 --interaction=1 --robust
[user@server ~/data_ex]$ ls
ex1_add.out.txt  ex2_domin.out.txt       ex3_add.out.txt  ex4_domin.out.txt       ex5_add.out.txt  ex6_domin.out.txt       ex.mldose  pheno.txt
ex2_2df.out.txt  ex2_over_domin.out.txt  ex4_2df.out.txt  ex4_over_domin.out.txt  ex6_2df.out.txt  ex6_over_domin.out.txt  ex.mlinfo
ex2_add.out.txt  ex2_recess.out.txt      ex4_add.out.txt  ex4_recess.out.txt      ex6_add.out.txt  ex6_recess.out.txt      ex.mlprob
