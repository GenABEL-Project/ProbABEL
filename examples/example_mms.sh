echo "analysis using MMSCORE"
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.mldose --sep="," -o mmscore --mmscore mmscore_InvSigma_aj.sex.age.dat
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.dose.fvi --sep="," -o mmscore_fv  --mmscore mmscore_InvSigma_aj.sex.age.dat
diff mmscore_add.out.txt mmscore_fv_add.out.txt
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.mlprob --ngpreds=2 --sep="," -o mmscore_prob --mmscore mmscore_InvSigma_aj.sex.age.dat
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.prob.fvi --ngpreds=2 --sep="," -o mmscore_prob_fv  --mmscore mmscore_InvSigma_aj.sex.age.dat

