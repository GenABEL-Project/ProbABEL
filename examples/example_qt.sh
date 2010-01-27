echo "analysing QT"
echo "base analysis"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o height_base
../bin/palinear -p height.txt -d test.mldose_fvf.fvi -i test.mlinfo -m test.map -c 19 -o height_base_fvf
echo "Option --allcov"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --allcov -o height_allcov
echo "Option --interaction=1"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --interaction=1 -o height_int1
echo "Option --robust"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --robust -o height_robust
echo "Option --robust --interaction=1"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --robust --interaction=1 -o height_robust_int1

echo "Option --ngp=2, mlprob file"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 -o height_ngp2
echo "Option --ngp=2 --allcov"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --allcov -o height_ngp2_allcov
echo "Option --ngp=2 --interaction=1"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --interaction=1 -o height_ngp2_int1
echo "Option --ngp=2 --robust"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --robust -o height_ngp2_robust
echo "Option --ngp=2 --robust --interaction=1"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --robust --interaction=1 -o height_ngp2_robust_int1

echo "Test with missing data"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 -o impute_with_missing
