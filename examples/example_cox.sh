echo "analysing Cox model"
../bin/pacoxph -p coxph_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o coxph
../bin/pacoxph -p coxph_data.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 -o coxph_fv
diff coxph_add.out.txt coxph_fv_add.out.txt
