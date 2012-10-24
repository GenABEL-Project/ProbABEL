echo "analysing Cox model"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph

# ../src/pacoxph \
#     -p ${srcdir}/coxph_data.txt \
#     -d ${srcdir}/test.dose.fvi \
#     -i ${srcdir}/test.mlinfo \
#     -m ${srcdir}/test.map \
#     -c 19 \
#     -o coxph_fv

# echo "pacoxph check: dose vs. dose_fv"
# diff coxph_add.out.txt coxph_fv_add.out.txt


# ../src/pacoxph \
#     -p ${srcdir}/coxph_data.txt \
#     -d ${srcdir}/test.mlprob \
#     -i ${srcdir}/test.mlinfo \
#     -m ${srcdir}/test.map \
#     --ngpreds=2 \
#     -c 19 \
#     -o coxph_prob

# echo "pacoxph check: dose vs. prob"
# diff coxph_add.out.txt coxph_prob_add.out.txt
