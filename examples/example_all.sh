

echo "extracting IDS..."
perl ../bin/extIDS.pl < test.mldose > mldose.IDS

echo "preparing phenofile"
R --no-save < ../bin/prepare_data.R

echo "analysing QT"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o height.out.txt

echo "analysing BT"
../bin/palogist -p logist_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o logist.out.txt

echo "analysing Cox model"
../bin/pacoxph -p coxph_data.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o coxph.out.txt
