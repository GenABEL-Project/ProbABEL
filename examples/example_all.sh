echo "extracting IDS..."
perl ../bin/extIDS.pl < test.mldose > mldose.IDS

echo "preparing phenofile"
R --no-save < ../bin/prepare_data.R

sh example_qt.sh
sh example_bt.sh
sh example_cox.sh
sh example_mms.sh
