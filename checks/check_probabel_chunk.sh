#!/bin/bash
# L.C. Karssen
# This script is used to test whether the probabel script works
# correctly when input is cut up in chunks

echo "Testing probabel..."

# Exit with error when one of the steps in the script fails
set -e

# -------- Set some default paths and file names -------
if [ -z ${srcdir} ]; then
    srcdir="."
fi
inputdir="${srcdir}/inputfiles/"
padir="${srcdir}/../src/"
results="${srcdir}/verified_results/"

dosefile="$inputdir/test.mldose"
probfile="$inputdir/test.mlprob"
infofile="$inputdir/test.mlinfo"
mapfile="$inputdir/test.map"
phenofile="$inputdir/height.txt"

probabel="${padir}/probabel"
probabelcfg="${padir}/probabel_config.cfg.example"
chunksep="_._chunk_._"
chrsep="_._chr_._"

# ------ Prepare probabel and the config file ------
sed 's;"./";"../src/";g' $probabel > probabel
chmod a+x probabel
cp $probabelcfg probabel_config.cfg
chmod +w probabel_config.cfg # Need this for make distcheck to work
cp $phenofile height.PHE

base="chr${chrsep}"
echo "TestCohortNoChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg

base="chr${chrsep}.chunk${chunksep}"
echo "TestCohortChunk,$base.info,$base.dose,$base.prob,$base.map" \
    >> probabel_config.cfg


# ---------- function definitions ----------
prepare_input ()
{
    if [ "$1" = "nochunk" ]; then
        # ------------------ No chunks test -------------------
        # Split the dose, prob and info files up into two chromosomes
        # with some  chunks
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

        WithOrWithout="without"
    elif [ "$1" = "chunk" ]; then
        # ------------------ Chunks test ----------------------
        # Split the dose and info files up into two chromosomes with
        # some chunks
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

        WithOrWithout="with"
    else
        echo "Run this function with one of these arguments: 'chunk'
        or 'nochunk'."
        exit 1
    fi

}


run_test ()
{
    # Run an analysis on dosage data
    outfile="height_add.out.txt"

    echo "Checking output using dosages $WithOrWithout chunks..."
    ./probabel 1 2 linear $1 --additive height > /dev/null

    blanks="                                                          "
    echo -n "  Verifying "
    if diff $outfile $results/$outfile; then
        echo -e "${outfile}${blanks:${#outfile}} OK"
    else
        echo -e "${outfile}${blanks:${#outfile}} FAILED"
        exit 1
    fi

    # Run an analysis on probabilities
    outfilelist="height_ngp2_2df.out.txt height_ngp2_recess.out.txt
    height_ngp2_over_domin.out.txt height_ngp2_domin.out.txt
    height_ngp2_add.out.txt"

    echo "Checking output using probabilities $WithOrWithout chunks..."
    ./probabel 1 2 linear $1 --allmodels height -o height_ngp2 > /dev/null
    for file in $outfilelist; do
        echo -n "  Verifying "
        if diff $file $results/$file; then
            echo -e "${file}${blanks:${#file}} OK"
        else
            echo -e "${file}${blanks:${#file}} FAILED"
            exit 1
        fi
    done
}

# ---------- Continuation of the main function ----------
prepare_input nochunk

run_test TestCohortNoChunk
echo "-------------------- Finished check without chunks --------------------"

prepare_input chunk

run_test TestCohortChunk
echo "-------------------- Finished check with chunks --------------------"
