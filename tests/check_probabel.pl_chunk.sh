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
results="${srcdir}/known_good_results/"

dosefile="$exampledir/test.mldose"
infofile="$exampledir/test.mlinfo"
mapfile="$exampledir/test.map"
phenofile="$exampledir/height.txt"
outfile="height_add.out.txt"

probabel="${padir}/probabel.pl"
probabelcfg="${padir}/probabel_config.cfg.example"
chunksep="_._chunk_._"
chrsep="_._chr_._"

# Prepare probabel.pl and the config file
sed 's;"./";"../src/";g' $probabel > probabel.pl
chmod a+x probabel.pl
cp $probabelcfg probabel_config.cfg
chmod +w probabel_config.cfg # Need this for make distcheck to work
cp $phenofile height.PHE

base="chr${chrsep}"
echo "TestCohortNoChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg

base="chunk${chunksep}.chr${chrsep}"
echo "TestCohortChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg


# ------------------ No chunks test -------------------
rm -f $outfile

# Split the dose and info files up into two chromosomes with some
# chunks
cut -d" " -f1,2,3,4 $dosefile > chr1.dose
cut -d" " -f1,2,5-7 $dosefile > chr2.dose

sed -n '1,3p' $infofile >  chr1.info
sed -n '1p'   $infofile >  chr2.info
sed -n '4,6p' $infofile >> chr2.info

sed -n '1,3p' $mapfile > chr1.map
sed -n '1p'   $mapfile >  chr2.map
sed -n '4,6p' $mapfile >> chr2.map

# Run an analysis
./probabel.pl 1 2 linear TestCohortNoChunk --additive height

# Final check:
echo "Checking output without chunks:"
diff $outfile $results/$outfile


# ------------------ Chunks test ----------------------
rm -f $outfile

# Split the dose and info files up into two chromosomes with some
# chunks
cut -d" " -f1,2,3   $dosefile > chunk1.chr1.dose
cut -d" " -f1,2,4   $dosefile > chunk2.chr1.dose
cut -d" " -f1,2,5,6 $dosefile > chunk1.chr2.dose
cut -d" " -f1,2,7   $dosefile > chunk2.chr2.dose

sed -n '1,2p' $infofile >  chunk1.chr1.info
sed -n '1p'   $infofile >  chunk2.chr1.info
sed -n '3p'   $infofile >> chunk2.chr1.info
sed -n '1p'   $infofile >  chunk1.chr2.info
sed -n '4,5p' $infofile >> chunk1.chr2.info
sed -n '1p'   $infofile >  chunk2.chr2.info
sed -n '6p'   $infofile >> chunk2.chr2.info

sed -n '1,2p' $mapfile >  chunk1.chr1.map
sed -n '1p'   $mapfile >  chunk2.chr1.map
sed -n '3p'   $mapfile >> chunk2.chr1.map
sed -n '1p'   $mapfile >  chunk1.chr2.map
sed -n '4,5p' $mapfile >> chunk1.chr2.map
sed -n '1p'   $mapfile >  chunk2.chr2.map
sed -n '6p'   $mapfile >> chunk2.chr2.map

# Run an analysis
./probabel.pl 1 2 linear TestCohortChunk --additive height

# Final check:
echo "Checking output with chunks:"
diff $outfile $results/$outfile
