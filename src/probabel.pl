#! /usr/bin/perl

#============================================================================
#
#           Filename:  probabel.pl
#
#        Description: Handy perl wrapper for ProbABEL functions
#
#            Version:  1.4
#            Created:  28-Oct-2008
#           Revision:  none
#  last modification:  07-Jan-2011
#
#             Author:  Lennart C. Karssen,
#                      Maksim V. Struchalin,
#                      Yurii S. Aulchenko
#            Company:  Erasmus MC,
#                      Epidemiology & Biostatistics Department,
#                      The Netherlands.
#              Email:  l.karssen@erasmusmc.nl,
#                      m.struchalin@erasmusmc.nl,
#                      i.aoultchenko@erasmusmc.nl
#
#============================================================================




$version="1.4";
$release_data="07-Jan-2011";

$_2df_file_postfix = "_2df.out.txt";
$_add_file_postfix = "_add.out.txt";
$_domin_file_postfix = "_domin.out.txt";
$_recess_file_postfix = "_recess.out.txt";
$_over_domin_file_postfix = "_over_domin.out.txt";




$separator_cfg = ",";
$separator_filename = "_._chr_._";

$base_path = "./";
$config = $base_path . "probabel_config.cfg";

@method = ("linear","logistic","coxph");
@anprog = ($base_path . "palinear", $base_path . "palogist", $base_path . "pacoxph");

%cohorts;
@mlinfo;
@mldose;
@mlprobe;
@legend;



#read config file
open(CFG,"$config") or die "$!";

<CFG>; #skip the first line (header)

for(my $i=0 ; my $line = <CFG> ; $i++)
	{
	next if (/^#/);
	chomp($line);
	@line_array = split(/$separator_cfg/, $line);
#	foreach (@line_array) {print $_."\n"}
#	print "\n";
	$cohorts{$line_array[0]} = $i;
	$mlinfo[$i]  = $line_array[1];
	$mldose[$i]  = $line_array[2];
	$mlprobe[$i] = $line_array[3];
	$legend[$i]  = $line_array[4];
#	print "legend[$i]=".$legend[$i]."\n";
#	print "mlprobe[$i]=".$mlprobe[$i]."\n";
	}
close(CFG);


if(@ARGV<6 || $ARGV[0] eq "--help")
	{
	print "\nUsage:
	probabel.pl chrom-start chrom-end method cohort <--allmodels OR --additive> trait <other available keys of ProbABEL functions>\n";
	print "\n	* chrom-start - the first chromosome number, chrom-end - the last one";
	print "\n	* method can be ";
	foreach $me(@method) {print "\"".$me."\", "};
	print "\n	* use --allmodels if you need dominant, recessive and heterozygous models
	  and --additive if additive only\n";
	print "	* Available cohorts is ";
 	foreach $coh(keys %cohorts) {print "\"".$coh."\", "};
	print "\n	* example:
	  probabel.pl 1 22 linear \"ERF\" --additive filename1 filename2
	  (filename has to be saved as filename.PHE)\n\n";


	print "\n Details:\n";
	print "   The probabel.pl script is used for analysis of imputed data.
   First you have to create file with the phenotypes that you are going to use.
   The first column contains ids in special order, the second: the trait which you are going analyze, the others containd covariates.
   Like this =>

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

   Then run the following on the Linux command line:
   probabel.pl 1 22 \"method\" \"cohort\" --model name_of_file
   Change \"method\" \"cohort\" --model to appropriate values\n";
	print "\n	Version: $version";
	print "\n	Release data: $release_data";
	print "\n\n	Authors: Lennart Karssen   - l.karssen\@erasmusmc.nl,
 		 Maksim Struchalin - m.struchalin\@erasmusmc.nl,
		 Yurii Aulchenko   - i.aoultchenko\@erasmusmc.nl.\n\n";
	exit;
	}


$startchr = $ARGV[0];
$endchr = $ARGV[1];
$method = $ARGV[2];
$chohort = $ARGV[3];
$model = $ARGV[4];



die "error: chrom-start is > 22" if($startchr > 22 && $startchr != "X") ;
die "error: chrom-end is > 22" if($endchr > 22 && $endchr != "X");
die "error: chrom-start > chrom-end" if($startchr > $endchr);


$cohort_position = $cohorts{$chohort};


if(!defined($cohort_position))
{
print "\nerror: Wrong cohort name, \"$chohort\" is not an available cohort.
Available cohorts are ";
foreach $coh(keys %cohorts) {print "\"".$coh."\", "};
print "\n\n";
exit;
}


@mlinfo_split = split(/$separator_filename/, $mlinfo[$cohort_position]);
@mldose_split = split(/$separator_filename/, $mldose[$cohort_position]);
@mlprobe_split = split(/$separator_filename/, $mlprobe[$cohort_position]);
@legend_split = split(/$separator_filename/, $legend[$cohort_position]);



$passed = 0;
for (my $i=0;$i<@method;$i++) {
		if ($ARGV[2] eq $method[$i]) {
				$passed = 1;
				$prog = $anprog[$i];
		}
}
die "error: Wrong method. method has to be one of: @method\n" if (!$passed);





$phename = $ARGV[5];
$keys="";


for ($i=6;$i<@ARGV;$i++)
	{
	$keys = $keys.$ARGV[$i]." ";
	}
chop($keys);



$model_option_num=0;

if($model eq "--additive")
	{
	$mldose_probe_split1 = $mldose_split[0];
	$mldose_probe_split2 = $mldose_split[1];
	$model_option_num=1;
	}
elsif($model eq "--allmodels")
	{
	$mldose_probe_split1 = $mlprobe_split[0];
	$mldose_probe_split2 = $mlprobe_split[1];
	$model_option_num=2;
	}
else
	{
	die "error: Wrong key for model. You can use \"--additive\" or \"--allmodels\" only\n";
	}



print "Start...\n";

$chr = $startchr;

# Separate command for the X chromosome.
if ($chr eq "X") {
    print `$prog -p $phename.PHE --ngpreds $model_option_num -i $mlinfo_split[0]$mlinfo_split[1] -d $mldose_probe_split1$mldose_probe_split2 -m $legend_split[0]$legend_split[1] --chrom $chr -o $phename $keys`;
    exit;
}

# Commands for the autosomes
print `$prog -p $phename.PHE --ngpreds $model_option_num -i $mlinfo_split[0]${chr}$mlinfo_split[1] -d $mldose_probe_split1${chr}$mldose_probe_split2 -m $legend_split[0]${chr}$legend_split[1] --chrom $chr -o $phename $keys`;
	for($chr=($startchr+1);$chr<=$endchr;$chr++)
		{
		print `$prog -p $phename.PHE --ngpreds $model_option_num -i $mlinfo_split[0]${chr}$mlinfo_split[1] -d $mldose_probe_split1${chr}$mldose_probe_split2 -m $legend_split[0]${chr}$legend_split[1] --chrom $chr -o $phename.$chr --no-head $keys`;

		if($model_option_num==2)
			{
			`cat $phename.${chr}$_2df_file_postfix >> ${phename}${_2df_file_postfix}`;
			`rm $phename.${chr}$_2df_file_postfix`;

			`cat $phename.${chr}$_add_file_postfix >> ${phename}${_add_file_postfix}`;
			`rm $phename.${chr}$_add_file_postfix`;

			`cat $phename.${chr}$_domin_file_postfix >> ${phename}${_domin_file_postfix}`;
			`rm $phename.${chr}$_domin_file_postfix`;

			`cat $phename.${chr}$_recess_file_postfix >> ${phename}${_recess_file_postfix}`;
			`rm $phename.${chr}$_recess_file_postfix`;

			`cat $phename.${chr}$_over_domin_file_postfix >> ${phename}${_over_domin_file_postfix}`;
			`rm $phename.${chr}$_over_domin_file_postfix`;
			}
		else
			{
			`cat $phename.${chr}$_add_file_postfix >> $phename${_add_file_postfix}`;
			print "cat $phename.${chr}$_add_file_postfix >> $phename${_add_file_postfix}\n";
			`rm $phename.${chr}$_add_file_postfix`;
			print "rm $phename.${chr}$_add_file_postfix\n";
			}
		}

print "Finished all chromosomes ...\n";
