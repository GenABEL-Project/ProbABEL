echo "analysing Cox model"

../bin/pacoxph -p coxph_data.txt -d test.mldose.2 -i test.mlinfo -m test.map -c 19 -o coxph.out.txt.dose

../bin/pacoxph -p coxph_data.txt -d test.mlprob --ngpreds=2 -i test.mlinfo -m test.map -c 19 -o coxph.out.txt.prob


diff coxph.out.txt.dose_add.out.txt coxph.out.txt.prob_add.out.txt
