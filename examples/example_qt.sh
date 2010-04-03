echo "analysing QT"

echo "base analysis"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 -o height_base
../bin/palinear -p height.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 -o height_base_fv
diff height_base_add.out.txt height_base_fv_add.out.txt

echo "Option --allcov"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --allcov -o height_allcov
../bin/palinear -p height.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 --allcov -o height_allcov_fv
diff height_allcov_add.out.txt height_allcov_fv_add.out.txt

echo "Option --interaction=1"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --interaction=1 -o height_int1
../bin/palinear -p height.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 --interaction=1 -o height_int1_fv
diff height_int1_add.out.txt height_int1_fv_add.out.txt

echo "Option --robust"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --robust -o height_robust
../bin/palinear -p height.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 --robust -o height_robust_fv
diff height_robust_add.out.txt height_robust_fv_add.out.txt

echo "Option --robust --interaction=1"
../bin/palinear -p height.txt -d test.mldose -i test.mlinfo -m test.map -c 19 --robust --interaction=1 -o height_robust_int1
../bin/palinear -p height.txt -d test.dose.fvi -i test.mlinfo -m test.map -c 19 --robust --interaction=1 -o height_robust_int1_fv
diff height_robust_int1_add.out.txt height_robust_int1_fv_add.out.txt

echo "Option --ngp=2, mlprob file"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 -o height_ngp2
../bin/palinear -p height.txt -d test.prob.fvi -i test.mlinfo -m test.map -c 19 --ngp=2 -o height_ngp2_fv
diff height_ngp2_add.out.txt height_ngp2_fv_add.out.txt
diff height_ngp2_domin.out.txt height_ngp2_fv_domin.out.txt
diff height_ngp2_over_domin.out.txt height_ngp2_fv_over_domin.out.txt
diff height_ngp2_recess.out.txt height_ngp2_fv_recess.out.txt
diff height_ngp2_2df.out.txt height_ngp2_fv_2df.out.txt

echo "Option --ngp=2 --allcov"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --allcov -o height_ngp2_allcov
../bin/palinear -p height.txt -d test.prob.fvi -i test.mlinfo -m test.map -c 19 --ngp=2 --allcov -o height_ngp2_allcov_fv
diff height_ngp2_allcov_add.out.txt height_ngp2_allcov_fv_add.out.txt
diff height_ngp2_allcov_domin.out.txt height_ngp2_allcov_fv_domin.out.txt
diff height_ngp2_allcov_over_domin.out.txt height_ngp2_allcov_fv_over_domin.out.txt
diff height_ngp2_allcov_recess.out.txt height_ngp2_allcov_fv_recess.out.txt
diff height_ngp2_allcov_2df.out.txt height_ngp2_allcov_fv_2df.out.txt

echo "Option --ngp=2 --interaction=1"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --interaction=1 -o height_ngp2_int1
../bin/palinear -p height.txt -d test.prob.fvi -i test.mlinfo -m test.map -c 19 --ngp=2 --interaction=1 -o height_ngp2_int1_fv
diff height_ngp2_int1_add.out.txt height_ngp2_int1_fv_add.out.txt
diff height_ngp2_int1_domin.out.txt height_ngp2_int1_fv_domin.out.txt
diff height_ngp2_int1_over_domin.out.txt height_ngp2_int1_fv_over_domin.out.txt
diff height_ngp2_int1_recess.out.txt height_ngp2_int1_fv_recess.out.txt
diff height_ngp2_int1_2df.out.txt height_ngp2_int1_fv_2df.out.txt

echo "Option --ngp=2 --robust"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --robust -o height_ngp2_robust
../bin/palinear -p height.txt -d test.prob.fvi -i test.mlinfo -m test.map -c 19 --ngp=2 --robust -o height_ngp2_robust_fv
diff height_ngp2_robust_add.out.txt height_ngp2_robust_fv_add.out.txt
diff height_ngp2_robust_domin.out.txt height_ngp2_robust_fv_domin.out.txt
diff height_ngp2_robust_over_domin.out.txt height_ngp2_robust_fv_over_domin.out.txt
diff height_ngp2_robust_recess.out.txt height_ngp2_robust_fv_recess.out.txt
diff height_ngp2_robust_2df.out.txt height_ngp2_robust_fv_2df.out.txt

echo "Option --ngp=2 --robust --interaction=1"
../bin/palinear -p height.txt -d test.mlprob -i test.mlinfo -m test.map -c 19 --ngp=2 --robust --interaction=1 -o height_ngp2_robust_int1
../bin/palinear -p height.txt -d test.prob.fvi -i test.mlinfo -m test.map -c 19 --ngp=2 --robust --interaction=1 -o height_ngp2_robust_int1_fv
diff height_ngp2_robust_int1_add.out.txt height_ngp2_robust_int1_fv_add.out.txt
diff height_ngp2_robust_int1_domin.out.txt height_ngp2_robust_int1_fv_domin.out.txt
diff height_ngp2_robust_int1_over_domin.out.txt height_ngp2_robust_int1_fv_over_domin.out.txt
diff height_ngp2_robust_int1_recess.out.txt height_ngp2_robust_int1_fv_recess.out.txt
diff height_ngp2_robust_int1_2df.out.txt height_ngp2_robust_int1_fv_2df.out.txt

