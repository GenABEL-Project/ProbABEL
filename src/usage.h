void print_usage (char * program_name, int exit_code)
{
	fprintf(stdout,"Usage: %s options\n",program_name);
	fprintf(stdout,"Options:\n");
	fprintf(stdout,"\t --pheno   : phenotype file name\n");
	fprintf(stdout,"\t --info    : information (e.g. MLINFO) file name\n");
	fprintf(stdout,"\t --dose    : predictor (e.g. MLDOSE/MLPROB) file name\n");
	fprintf(stdout,"\t --map     : [optional] map file name\n");
	fprintf(stdout,"\t --nids    : [optional] number of people to analyse\n");
	fprintf(stdout,"\t --chrom   : [optional] chromosome (to be passed to output)\n");
	fprintf(stdout,"\t --out     : [optional] output file name (default is regression.out.txt)\n");
	fprintf(stdout,"\t --skipd   : [optional] how many columns to skip in the predictor\n\t              (dose/prob) file (default 2)\n");
#if COXPH
	fprintf(stdout,"\t --ntraits : [optional] how many traits are analysed (default 2)\n");
#else
	fprintf(stdout,"\t --ntraits : [optional] how many traits are analysed (default 1)\n");
#endif
	fprintf(stdout,"\t --ngpreds : [optional] how many predictor columns per marker\n\t              (default 1 = MLDOSE; else use 2 for MLPROB)\n");
	fprintf(stdout,"\t --separat : [optional] character to separate fields (default is space)\n");
	fprintf(stdout,"\t --score   : use score test\n");
	fprintf(stdout,"\t --no-head : do not report header line\n");
	fprintf(stdout,"\t --allcov  : report estimates for all covariates (large outputs!)\n");
	fprintf(stdout,"\t --interaction: Which covariate to use for interaction with SNP analysis (default is no interaction, 0)\n");
	fprintf(stdout,"\t --interaction_only: like previous but without covariate acting in interaction with SNP (default is no interaction, 0)\n");
	fprintf(stdout,"\t --mmscore : score test in samples of related individuals. File with inverse of variance-covariance matrix (for palinear) or inverse correlation (for palogist) as input parameter\n");
	fprintf(stdout,"\t --robust  : report robust (aka sandwich, aka Hubert-White) standard errors\n");
	fprintf(stdout,"\t --help    : print help\n");
	exit(exit_code);
}

void print_help (char * program_name, int exit_code)
{
	print_usage(program_name,exit_code);
}
