echo "analysing Linear model"

../bin/palinear -p height.txt -d test.mldose.2 -i test.mlinfo -m test.map -c 19 -o height.out.txt.dose

../bin/palinear -p height.txt -d test.mlprob --ngpreds=2 -i test.mlinfo -m test.map -c 19 -o height.out.txt.prob

diff height.out.txt.dose_add.out.txt height.out.txt.prob_add.out.txt
