#!/bin/bash
# L.C. Karssen
# This script is used to test whether probabel.pl works correctly when
# input is cut up in chunks

echo "Testing probabel.pl..."

if [ -z ${srcdir} ]; then
    srcdir="."
fi
exampledir="${srcdir}/../examples/"
padir="${srcdir}/../src/"
results="${srcdir}/verified_results/"

dosefile="$exampledir/test.mldose"
probfile="$exampledir/test.mlprob"
infofile="$exampledir/test.mlinfo"
mapfile="$exampledir/test.map"
phenofile="$exampledir/height.txt"

probabel="${padir}/probabel.pl"
probabelcfg="${padir}/probabel_config.cfg.example"
chunksep="_._chunk_._"
chrsep="_._chr_._"

# ------ Prepare probabel.pl and the config file ------
sed 's;"./";"../src/";g' $probabel > probabel.pl
chmod a+x probabel.pl
cp $probabelcfg probabel_config.cfg
chmod +w probabel_config.cfg # Need this for make distcheck to work
cp $phenofile height.PHE

base="chr${chrsep}"
echo "TestCohortNoChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg

base="chr${chrsep}.chunk${chunksep}"
echo "TestCohortChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg


# ------------------ No chunks test -------------------
# Split the dose, prob and info files up into two chromosomes with some
# chunks
awk '{print $1,$2,$3,$4}'    $dosefile > chr1.dose
awk '{print $1,$2,$5,$6,$7}' $dosefile > chr2.dose

awk '{print $1,$2,$3,$4,$5,$6}'          $probfile > chr1.prob
awk '{print $1,$2,$7,$8,$9,$10,$11,$12}' $probfile > chr2.prob

sed -n '1,3p' $infofile >  chr1.info
sed -n '1p'   $infofile >  chr2.info
sed -n '4,6p' $infofile >> chr2.info

sed -n '1,3p' $mapfile > chr1.map
sed -n '1p'   $mapfile >  chr2.map
sed -n '4,6p' $mapfile >> chr2.map

# Run an analysis on dosage data
outfile="height_add.out.txt"
rm -f $outfile
./probabel.pl 1 2 linear TestCohortNoChunk --additive height
echo "Checking output using dosages without chunks..."
diff $outfile $results/$outfile


# Run an analysis on probabilities
outfile="height_2df.out.txt"
rm -f $outfile
./probabel.pl 1 2 linear TestCohortNoChunk --allmodels height
echo "Checking output using probabilities without chunks..."
diff $outfile $results/$outfile



# ------------------ Chunks test ----------------------
# Split the dose and info files up into two chromosomes with some
# chunks
awk '{print $1,$2,$3}'    $dosefile > chr1.chunk1.dose
awk '{print $1,$2,$4}'    $dosefile > chr1.chunk2.dose
awk '{print $1,$2,$5,$6}' $dosefile > chr2.chunk1.dose
awk '{print $1,$2,$7}'    $dosefile > chr2.chunk2.dose

awk '{print $1,$2,$3,$4}'        $probfile > chr1.chunk1.prob
awk '{print $1,$2,$5,$6}'        $probfile > chr1.chunk2.prob
awk '{print $1,$2,$7,$8,$9,$10}' $probfile > chr2.chunk1.prob
awk '{print $1,$2,$11,$12}'      $probfile > chr2.chunk2.prob

sed -n '1,2p' $infofile >  chr1.chunk1.info
sed -n '1p'   $infofile >  chr1.chunk2.info
sed -n '3p'   $infofile >> chr1.chunk2.info
sed -n '1p'   $infofile >  chr2.chunk1.info
sed -n '4,5p' $infofile >> chr2.chunk1.info
sed -n '1p'   $infofile >  chr2.chunk2.info
sed -n '6p'   $infofile >> chr2.chunk2.info

sed -n '1,2p' $mapfile >  chr1.chunk1.map
sed -n '1p'   $mapfile >  chr1.chunk2.map
sed -n '3p'   $mapfile >> chr1.chunk2.map
sed -n '1p'   $mapfile >  chr2.chunk1.map
sed -n '4,5p' $mapfile >> chr2.chunk1.map
sed -n '1p'   $mapfile >  chr2.chunk2.map
sed -n '6p'   $mapfile >> chr2.chunk2.map

# Run an analysis on dosage data
outfile="height_add.out.txt"
rm -f $outfile
./probabel.pl 1 2 linear TestCohortChunk --additive height
echo "Checking output using dosages with chunks..."
diff $outfile $results/$outfile


# Run an analysis on probabilities
outfile="height_2df.out.txt"
rm -f $outfile
./probabel.pl 1 2 linear TestCohortNoChunk --allmodels height
echo "Checking output using probabilities without chunks:"
diff $outfile $results/$outfile
