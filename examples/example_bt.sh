echo "analysing BT"
../bin/palogist -p logist_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o logist
../bin/palogist -p logist_data.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 -o logist_fv
diff logist_add.out.txt logist_fv_add.out.txt

