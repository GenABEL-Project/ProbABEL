echo "analysing Cox model"
../bin/pacoxph -p coxph_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o coxph.out.txt
