echo "analysing Cox model"
../bin/pacoxph -p coxph_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o coxph.out.txt
../bin/pacoxph -p coxph_data.txt -d test.mldose_fvf.fvi -i test.mlinfo -m test.map -c 19 -o coxph.out.txt_fvf
