echo "analysing BT"
../bin/palogist -p logist_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o logist.out.txt
../bin/palogist -p logist_data.txt -d test.mldose_fvf.fvi -i test.mlinfo -m test.map -c 19 -o logist.out.txt_fvf
