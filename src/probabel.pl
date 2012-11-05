#! /usr/bin/perl -W
#==========================================================================
#
#           Filename:  probabel.pl
#
#        Description: Handy perl wrapper for ProbABEL functions
#
#==========================================================================
use strict;

#==========================================================================
# Set variables
my $version="PROBABEL_VERSION";

# Define some filename postfixes
my $_2df_file_postfix = "_2df.out.txt";
my $_add_file_postfix = "_add.out.txt";
my $_domin_file_postfix = "_domin.out.txt";
my $_recess_file_postfix = "_recess.out.txt";
my $_over_domin_file_postfix = "_over_domin.out.txt";

# Separators in the config file
my $separator_cfg = ",";
my $chr_replacement = "_._chr_._";
my $chunk_replacement = "_._chunk_._";

# Set file locations
my $base_path = "./";
my @anprog = ($base_path . "palinear",
	      $base_path . "palogist",
	      $base_path . "pacoxph");
my $config = "probabel_config.cfg";

# Define the regression methods that are implemented
my @method = ("linear", "logistic", "coxph");

my %cohorts;
my @mlinfos;
my @mldoses;
my @mlprobs;
my @legends;


#==========================================================================
# Read config file
open(CFG, "$config") or die "Reading configuration file $config failed: $!" .
    "\nDid you forget to edit and rename the probabel_config.cfg.example file?\n";

<CFG>; #skip the first line (header)

for(my $i=0 ; my $line = <CFG> ; $i++)
{
    chomp($line);
    next if ($line =~ /^#/);
    next if ($line =~ /^$/);
    my @line_array = split(/$separator_cfg/, $line);
    $cohorts{$line_array[0]} = $i;
    $mlinfos[$i]  = $line_array[1];
    $mldoses[$i]  = $line_array[2];
    $mlprobs[$i] = $line_array[3];
    $legends[$i]  = $line_array[4];
}
close(CFG);


#==========================================================================
# Print usage info if arguments are not correct
if(@ARGV<6 || $ARGV[0] eq "--help" || $ARGV[0] eq "-h") {
    print "Usage:
	probabel.pl chrom-start chrom-end method cohort <--allmodels OR --additive> trait <other available options of ProbABEL functions>\n";
    print "\n	* chrom-start - the first chromosome number, chrom-end - the last one; X or Y have to be run separately (specify them twice, once as chrom-start and once as chrom-end)";
    print "\n	* method can be ";
    foreach my $me(@method) {print "\"".$me."\", "};
    print "\n	* use --allmodels if you need dominant, recessive and heterozygous models
	  and --additive if additive only\n";
    print "	* Available cohorts are ";
    foreach my $coh(keys %cohorts) {
	print "\"".$coh."\", "
    };
    print "\n	* example:
	  probabel.pl 1 22 linear \"ERF\" --additive filename
	  (filename has to be saved as filename.PHE)\n\n";

    if(@ARGV == 1 && ($ARGV[0] eq "--help" || $ARGV[0] eq "-h")) {
	print "\nDetails:\n";
	print " The probabel.pl script is used for analysis of imputed data. First you have to create a file with the phenotype values that you are going to use. The first column contains ids in special order, the second one contains the trait which you are going	analyze, the others contain covariates.  For example:

	id         phen1 covariate1  covariate2
	1_2094  0     334         0
	1_5060  1     56          1
	1_4077  1     346         6
	.
	.
	.

	This implies the model:
	phen1 ~ covariate1 + covariate2 + SNP


	Then save it to folder where you are doing the analysis. The name of the file must be name_of_file.PHE, where name_of_file is any name.

	Then run the following on the command line:
	  probabel.pl 1 22 \"method\" \"cohort\" --model name_of_file
	Change \"method\", \"cohort\" and --model to appropriate values\n";
	print "\n	Version: $version";
	print "\n\n	Authors: Lennart Karssen   - l.karssen\@erasmusmc.nl,
 		 Maksim Struchalin - m.struchalin\@erasmusmc.nl,
		 Yurii Aulchenko   - yurii.aulchenko\@gmail.com.\n\n";
    }
    else {
	print "Type probabel.pl --help for more details.\n";
    }
    exit;
}


#==========================================================================
# Put the command line arguments into variables and verify them
my $startchr = $ARGV[0];
my $endchr = $ARGV[1];
my $method = $ARGV[2];
my $chohort = $ARGV[3];
my $model = $ARGV[4];

die "error: chrom-start is > 22" if($startchr > 22 && $startchr != "X") ;
die "error: chrom-end is > 22" if($endchr > 22 && $endchr != "X");
die "error: chrom-start > chrom-end" if($startchr > $endchr);


my $cohort_position = $cohorts{$chohort};


if(!defined($cohort_position)) {
    print "\nerror: Wrong cohort name, \"$chohort\" is not an available cohort.
Available cohorts are ";
    foreach my $coh(keys %cohorts) {
	print "\"".$coh."\", ";
    }
    print "\n\n";
    exit;
}

my $mlinfo = $mlinfos[$cohort_position];
my $mldose = $mldoses[$cohort_position];
my $mlprob = $mlprobs[$cohort_position];
my $legend = $legends[$cohort_position];


my $passed = 0;
my $prog;
for (my $i=0; $i<@method; $i++) {
    if ($ARGV[2] eq $method[$i]) {
	$passed = 1;
	$prog = $anprog[$i];
    }
}
die "error: Wrong method. method has to be one of: @method\n" if (!$passed);


my $phename = $ARGV[5];
my $outfile_prefix = $phename;
my $keys="";
for (my $i=6; $i<@ARGV; $i++) {
    if ($ARGV[$i] eq "-o")
    {
	# Apparently the user wants to change the output file name
	# Let's interpret this as an addition to our own prefix
	$outfile_prefix = $outfile_prefix.$ARGV[$i+1];

	# Skip the next argument (supposedly the addition to the
	# output file name).
	$i++;
    }
    else
    {
	$keys = $keys.$ARGV[$i]." ";
    }
}
chop($keys);

my $model_option_num = 0;
my $mldose_prob;
if($model eq "--additive") {
    $mldose_prob = $mldose;
    $model_option_num = 1;
} elsif($model eq "--allmodels") {
    $mldose_prob = $mlprob;
    $model_option_num = 2;
} else {
    die "error: Wrong key for model. You can use \"--additive\" or \"--allmodels\" only\n";
}


#==========================================================================
# Start the analysis now that the input has been validated
print "Start...\n";

my $chr = $startchr;
my $hadhead=0;
my $head;
my $mlinfo_arg;
my $mldose_arg;
my $legend_arg;

# Separate command for the sex chromosomes.
if ($chr eq "X" || $chr eq "Y") {
    $mlinfo_arg = $mlinfo;
    $mlinfo_arg =~ s/$chr_replacement/$chr/g;

    $mldose_arg = $mldose_prob;
    $mldose_arg =~ s/$chr_replacement/$chr/g;

    $legend_arg = $legend;
    $legend_arg =~ s/$chr_replacement/$chr/g;

    if($hadhead==0) {
	$head="";
	$hadhead=1;
    } else {
	my $head="--no-head";
    }

    system "$prog -p $phename.PHE --ngpreds $model_option_num -i $mlinfo_arg -d $mldose_arg -m $legend_arg --chrom $chr -o $outfile_prefix $head $keys";

    exit;
}

# Commands for the autosomes
for($chr=$startchr; $chr<=$endchr; $chr++) {

    my $nrchunks = 0;
    # Find out the number of chunks for the current chromosome
    my $infofiles = $mlinfo;
    $infofiles =~ s/$chr_replacement/$chr/g;
    $infofiles =~ s/$chunk_replacement/*/g;
    $nrchunks = `ls $infofiles 2>/dev/null | wc -l`;
    if ($nrchunks==0) {
	# If no chunked info files exist the 'wc -l' command returns 0
	# so that actually means 1 chunk containing all data.
	$nrchunks = 1;
    }
    print "Nr. of chunks: $nrchunks";

    # Loop over all chunks
    for (my $chunk=1; $chunk <= $nrchunks; $chunk++)
    {
	if($hadhead==0) {
	    $head="";
	    $hadhead=1;
	} else {
	    $head="--no-head";
	}
	$mlinfo_arg = $mlinfo;
	$mlinfo_arg =~ s/$chr_replacement/$chr/g;
	$mlinfo_arg =~ s/$chunk_replacement/$chunk/g;

	$mldose_arg = $mldose_prob;
	$mldose_arg =~ s/$chr_replacement/$chr/g;
	$mldose_arg =~ s/$chunk_replacement/$chunk/g;

	$legend_arg = $legend;
	$legend_arg =~ s/$chr_replacement/$chr/g;
	$legend_arg =~ s/$chunk_replacement/$chunk/g;

	my $command = "$prog -p $phename.PHE --ngpreds $model_option_num -i $mlinfo_arg -d $mldose_arg -m $legend_arg --chrom $chr -o $outfile_prefix.chunk$chunk.chr$chr $head $keys";
	print "$command \n";
	system $command;

	if($model_option_num==2)
	{
	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_2df_file_postfix >> ${outfile_prefix}.${chr}${_2df_file_postfix}`;
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_2df_file_postfix`;

	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix >> ${outfile_prefix}.${chr}${_add_file_postfix}`;
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix`;

	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_domin_file_postfix >> ${outfile_prefix}.${chr}${_domin_file_postfix}`;
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_domin_file_postfix`;

	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_recess_file_postfix >> ${outfile_prefix}.${chr}${_recess_file_postfix}`;
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_recess_file_postfix`;

	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_over_domin_file_postfix >> ${outfile_prefix}.${chr}${_over_domin_file_postfix}`;
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_over_domin_file_postfix`;
	} else {
	    `cat $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix >> $outfile_prefix.${chr}${_add_file_postfix}`;
	    print "cat $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix >> $outfile_prefix.chr${chr}${_add_file_postfix}\n";
	    `rm $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix`;
	    print "rm $outfile_prefix.chunk${chunk}.chr${chr}$_add_file_postfix\n";
	}
    }

    if($model_option_num==2)
    {
	`cat $outfile_prefix.${chr}$_2df_file_postfix >> ${outfile_prefix}${_2df_file_postfix}`;
	`rm $outfile_prefix.${chr}$_2df_file_postfix`;

	`cat $outfile_prefix.${chr}$_add_file_postfix >> ${outfile_prefix}${_add_file_postfix}`;
	`rm $outfile_prefix.${chr}$_add_file_postfix`;

	`cat $outfile_prefix.${chr}$_domin_file_postfix >> ${outfile_prefix}${_domin_file_postfix}`;
	`rm $outfile_prefix.${chr}$_domin_file_postfix`;

	`cat $outfile_prefix.${chr}$_recess_file_postfix >> ${outfile_prefix}${_recess_file_postfix}`;
	`rm $outfile_prefix.${chr}$_recess_file_postfix`;

	`cat $outfile_prefix.${chr}$_over_domin_file_postfix >> ${outfile_prefix}${_over_domin_file_postfix}`;
	`rm $outfile_prefix.${chr}$_over_domin_file_postfix`;
    } else {
	`cat $outfile_prefix.${chr}$_add_file_postfix >> $outfile_prefix${_add_file_postfix}`;
	print "cat $outfile_prefix.${chr}$_add_file_postfix >> $outfile_prefix${_add_file_postfix}\n";
	`rm $outfile_prefix.${chr}$_add_file_postfix`;
	print "rm $outfile_prefix.${chr}$_add_file_postfix\n";
    }
}

print "Finished all chromosomes ...\n";
