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
../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph_fv

diff coxph_add.out.txt coxph_fv_add.out.txt
