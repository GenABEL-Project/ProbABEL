echo "analysing BT"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/palogist -p ${srcdir}/logist_data.txt -d ${srcdir}/test.mldose -i ${srcdir}/test.mlinfo -m ${srcdir}/test.map -c 19 -o logist
../src/palogist -p ${srcdir}/logist_data.txt -d ${srcdir}/test.dose.fvi -i ${srcdir}/test.mlinfo -m ${srcdir}/test.map -c 19 -o logist_fv

diff logist_add.out.txt logist_fv_add.out.txt
