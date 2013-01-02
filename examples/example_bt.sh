echo "analysing BT"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o logist

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o logist_fv

echo "BT check: dose vs. dose_fv"
diff logist_add.out.txt logist_fv_add.out.txt

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob_fv

echo "BT check: prob vs. prob_fv"
diff logist_prob_add.out.txt logist_prob_fv_add.out.txt

# Commented out because of slightly different output formats. We need
# something smart here.
# echo "BT check: prob vs. dose"
# diff logist_prob_add.out.txt logist_add.out.txt
# echo "BT check: prob_fv vs. dose_fv"
# diff logist_prob_fv_add.out.txt logist_fv_add.out.txt
