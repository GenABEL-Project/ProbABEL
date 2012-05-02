echo "Analysing Cox model"

outfile_base="coxph.out.txt.dose"
outfile_dose="${outfile_base}_add.out.txt"
../bin/pacoxph -p coxph_data.txt -d test.mldose.2 -i test.mlinfo -m test.map -c 19 -o $outfile_base

outfile_base="coxph.out.txt.prob"
outfile_prob="${outfile_base}_add.out.txt"
../bin/pacoxph -p coxph_data.txt -d test.mlprob --ngpreds=2 -i test.mlinfo -m test.map -c 19 -o $outfile_base

# This won't work, because differences are in second/third decimal
#diff $outfile_dose $outfile_prob

tmp_dose=`tempfile`
tmp_prob=`tempfile`

gawk 'NR>1 {
       for (i=1; i<=8; i++)
       {
	  printf "%s ", $i
       };
       for (i=9; i<=NF; i++)
       {
	  printf "%6.2f ", $i
       };
       printf "\n"
    }' $outfile_dose > $tmp_dose

gawk 'NR>1 {
       for (i=1; i<=8; i++)
       {
	  printf "%s ", $i
       };
       for (i=9; i<=NF; i++)
       {
	  printf "%6.2f ", $i
       };
       printf "\n"
    }' $outfile_prob > $tmp_prob

diff $tmp_dose $tmp_prob
rm $tmp_dose $tmp_prob
