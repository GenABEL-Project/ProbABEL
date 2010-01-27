echo "analysis using MMSCORE"
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.mldose --sep="," -o mmscore.output  --mmscore mmscore_InvSigma_aj.sex.age.dat
../bin/palinear -p mmscore_pheno.PHE  -i mmscore_gen.mlinfo -d mmscore_gen.mldose_fvf.fvi --sep="," -o mmscore.output_fvf  --mmscore mmscore_InvSigma_aj.sex.age.dat
