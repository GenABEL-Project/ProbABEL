echo "analysing QT with missing"
echo "base analysis"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 -o impute_with_missing_base
echo "Option --allcov"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 --allcov -o impute_with_missing_allcov
echo "Option --interaction=1"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 --interaction=1 -o impute_with_missing_int1
echo "Option --robust"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 --robust -o impute_with_missing_robust
echo "Option --robust --interaction=1"
../bin/palinear -p impute_with_missing.PHE -d impute_with_missing.dose.fvi -i impute_with_missing.mlinfo -c 23 --robust --interaction=1 -o impute_with_missing_robust_int1

