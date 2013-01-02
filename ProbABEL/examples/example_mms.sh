echo "analysis using MMSCORE"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mldose \
    --sep="," \
    -o mmscore \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

echo "mmscore check: dose vs. dose_fv"
diff mmscore_add.out.txt mmscore_fv_add.out.txt


../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mlprob \
    --ngpreds=2 --sep="," \
    -o mmscore_prob \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.prob.fvi \
    --ngpreds=2 --sep="," \
    -o mmscore_prob_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

echo "mmscore check: prob vs. prob_fv"
diff mmscore_prob_add.out.txt mmscore_prob_fv_add.out.txt

# Commented out because of slightly different output formats. We need
# something smart here.
# echo "mmscore check: prob vs. dose"
# diff mmscore_prob_add.out.txt mmscore_add.out.txt

# echo "mmscore check: prob_fv vs. dose_fv"
# diff mmscore_prob_fv_add.out.txt mmscore_fv_add.out.txt
