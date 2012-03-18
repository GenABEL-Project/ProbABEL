#include <string>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include "mematrix.h"

#include "fvlib/AbstractMatrix.h"
#include "fvlib/CastUtils.h"
#include "fvlib/const.h"
#include "fvlib/convert_util.h"
#include "fvlib/FileVector.h"
#include "fvlib/frutil.h"
#include "fvlib/frversion.h"
#include "fvlib/Logger.h"
#include "fvlib/Transposer.h"

extern bool is_interaction_excluded;

void error(const char * format, ...)
{
    va_list args;
    char buffer[256];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    printf("ERROR: %s\n",buffer);
    exit(EXIT_FAILURE);
}

unsigned int Nmeasured(char * fname, int nphenocols, int npeople)
{
    int ncov = nphenocols - 2;
    int nids_all = npeople;

    FILE * infile;
    // first pass -- find unmeasured people
    if ((infile=fopen(fname,"r"))==NULL) {
	fprintf(stderr,"Nmeasured: cannot open file %s\n",fname);
    }
    char tmp[100];

    for (int i=0;i<nphenocols;i++)
    {
	fscanf(infile,"%s", tmp);
	//		printf("%s ",tmp);
    } 	//printf("\n");

    unsigned short int * allmeasured = new unsigned short int [npeople];
    int nids = 0;
    for (int i = 0;i<npeople;i++)
    {
	allmeasured[i] = 1;
	fscanf(infile,"%s", tmp);
	for (int j=1;j<nphenocols;j++)
	{
	    fscanf(infile,"%s", tmp);
	    if (tmp[0]=='N' || tmp[0]=='n') allmeasured[i]=0;
	}
	if (allmeasured[i]==1) nids++;
    }
    fclose(infile);
    return(nids);
}

class phedata
{
public:
    int nids_all;
    int nids;
    int noutcomes;
    int ncov;
    unsigned short int * allmeasured;
    mematrix<double> X;
    mematrix<double> Y;
    std::string * idnames;
    std::string model;
    std::string * model_terms;
    int n_model_terms;
    phedata(char * fname, int noutc, int npeople, int interaction, bool iscox)
    {
	static const unsigned int BFS = 1000;
	std::ifstream myfile(fname);
	char line[BFS];
	char tmp[100];
	noutcomes = noutc;

	int nphenocols=0;
	int savenpeople = npeople;
	npeople=0;
	if (myfile.is_open())
	{
	    myfile.getline(line,BFS);
	    std::stringstream line_stream(line);
	    //			std::cout << line << "\n ";
	    while (line_stream >> tmp)
	    {

		nphenocols++;
		//				std::cout << tmp << " " << nphenocols << " ";
	    }
	    while (myfile.getline(line,BFS)) {
		int tmplins = 0;
		std::stringstream line_stream(line);
		while (line_stream >> tmp)
		    tmplins++;
		if (tmplins != nphenocols)
		{
		    fprintf(stderr,"phenofile: number of variables different from %d in line %d\n",nphenocols,tmplins);
		    myfile.close();
		    exit(1);
		}
		npeople++;
	    };
	    myfile.close();
	}
	else
	{
	    fprintf(stderr,"Unable to open file %s\n",fname);
	    exit(1);
	}
	fprintf(stdout,"Actual number of people in phenofile = %d",npeople);
	if (savenpeople>0)
	{
	    npeople = savenpeople;
	    fprintf(stdout,"; using only %d first\n",npeople);
	}
	else
	{
	    fprintf(stdout,"; using all of these\n");
	}

	ncov = nphenocols - 1 - noutcomes;
	nids_all = npeople;
	model_terms = new std::string [ncov+2];

	FILE * infile;
	// first pass -- find unmeasured people
	if ((infile=fopen(fname,"r"))==NULL) {
	    fprintf(stderr,"phedata: cannot open file %s\n",fname);
	}


	fscanf(infile,"%s", tmp);
	model = "( ";
	fscanf(infile,"%s", tmp);
	model = model + tmp;
	for (int i = 1;i < noutcomes;i++)
	{
	    fscanf(infile,"%s", tmp);
	    model = model + " , ";
	    model = model + tmp;
	}
	n_model_terms = 0;
#if COXPH
	model = model + " ) ~ ";
#else
	model = model + " ) ~ mu + ";
	model_terms[n_model_terms++] = "mu";
#endif

	if (nphenocols>noutcomes+1)
	{
	    fscanf(infile,"%s", tmp);
	    model = model + tmp;
	    model_terms[n_model_terms++] = tmp;
	    for (int i=(2+noutcomes);i<nphenocols;i++)
	    {
		fscanf(infile,"%s", tmp);

		//				if(iscox && ) {if(n_model_terms+1 == interaction-1) {continue;} }
		//				else      {if(n_model_terms+1 == interaction) {continue;} }
		model = model + " + ";
		model = model + tmp;
		model_terms[n_model_terms++] = tmp;
	    }
	}
	model = model + " + SNP_A1";
	if(interaction!=0)
	{
	    if(iscox) {model = model + " + " +model_terms[interaction-1]+"*SNP_A1";}
	    else      {model = model + " + " +model_terms[interaction]+"*SNP_A1";}
	}
	model_terms[n_model_terms++] = "SNP_A1";

	if(is_interaction_excluded) // exclude covariates from covariate names
	{
	    if(iscox) {std::cout<<"model is running without "<<model_terms[interaction-1]<<", term\n";}
	    else      {std::cout<<"model is running without "<<model_terms[interaction]<<", term\n";}
	}


#if LOGISTIC
	std::cout << "Logistic ";
#elif LINEAR
	std::cout << "Linear ";
#elif COXPH
	std::cout << "Coxph ";
#else
	std::cout << "Unrecognised ";
#endif
	std::cout << "model: " << model << "\n";

	allmeasured = new unsigned short int [npeople];
	nids = 0;
	for (int i = 0;i<npeople;i++)
	{
	    allmeasured[i] = 1;
	    for (int j=0;j<nphenocols;j++)
	    {
		fscanf(infile,"%s", tmp);
		if (j>0 && (tmp[0]=='N' || tmp[0]=='n')) allmeasured[i]=0;
	    }
	    if (allmeasured[i]==1) nids++;
	}
	fclose(infile);
	//		printf("npeople = %d, no. all measured = %d\n",nids_all,nids);

	// allocate objects
	int ntmpcov = 1;
	if (ncov>0) ntmpcov = ncov;
	idnames = new std::string [nids];
	X.reinit(nids,ntmpcov);
	Y.reinit(nids,noutcomes);

	// second pass -- read the data
	if ((infile=fopen(fname,"r"))==NULL) {
	    fprintf(stderr,"phedata: cannot open file %s\n",fname);
	    exit(1);
	}

	for (int i=0;i<nphenocols;i++)
	{
	    fscanf(infile,"%s", tmp);
	}

	int k =0;
	int m =0;
	for (int i = 0;i<npeople;i++)
	    if (allmeasured[i]==1)
	    {
		fscanf(infile,"%s", tmp);
		idnames[m] = tmp;
		for (int j=0;j<noutcomes;j++)
		{
		    fscanf(infile,"%s", tmp);
		    Y.put(atof(tmp),m,j);
		}
		for (int j=(1+noutcomes);j<nphenocols;j++)
		{
		    fscanf(infile,"%s", tmp);
		    X.put(atof(tmp),m,(j-1-noutcomes));
		}
		m++;
	    }
	    else
		for (int j=0;j<nphenocols;j++) fscanf(infile,"%s", tmp);
	fclose(infile);
    }
    ~phedata()
    {
	//		delete X;
	//		delete Y;
	//		delete [] allmeasured;
    }
};

class gendata
{
public:
    int nsnps;
    int nids;
    int ngpreds;
    gendata();
    void re_gendata(char * fname, int insnps, int ingpreds,
		    int npeople, int nmeasured,
		    unsigned short int * allmeasured,
		    int skipd,
		    std::string * idnames);
    void re_gendata(string filename, int insnps, int ingpreds,
		    int npeople, int nmeasured,
		    unsigned short int * allmeasured,
		    std::string * idnames);
    ~gendata();
    void get_var(int var, float * data);
    // MAKE THAT PRIVATE, ACCESS THROUGH GET_SNP
    // ANOTHER PRIVATE OBJECT IS A POINTER TO DATABELBASECPP
    // UPDATE SNP, ALL REGRESSION METHODS: ACCOUNT FOR MISSING
private:
    mematrix<float> G;
    AbstractMatrix * DAG;
    unsigned short int * DAGmask;
    //	mematrix<double> G;
};

void gendata::get_var(int var, float * data)
{
    if (DAG == NULL)
	for (int i=0;i<G.nrow;i++) data[i] = G.get(i,var);
    else if (DAG != NULL) {
	float tmpdata[DAG->getNumObservations()];
	DAG->readVariableAs((unsigned long int) var, tmpdata);
	unsigned int j = 0;
	for (unsigned int i=0;i<DAG->getNumObservations();i++) if (!DAGmask[i]) data[j++] = tmpdata[i];
	//fprintf(stdout,"%i %i %i\n",j,DAG->get_nobservations(),nids);
    }
    else error("cannot get gendata");
}

gendata::gendata()
{
    nsnps=nids=ngpreds=0;
}

void gendata::re_gendata(string filename, int insnps, int ingpreds, int npeople,
			 int nmeasured, unsigned short int * allmeasured,
			 std::string * idnames)
{
    nsnps = insnps;
    ngpreds = ingpreds;
    DAG = new FileVector(filename,128,true);
    DAGmask = new unsigned short int [DAG->getNumObservations()];
    if (DAG->getNumObservations() != npeople) error("dimension of fvf-data and phenotype data do not match\n");
    if (DAG->getNumVariables() != insnps*ingpreds) error("dimension of fvf-data and mlinfo data do not match\n");
    long int j = -1;
    for (unsigned int i=0;i<npeople;i++)
    {
	if (allmeasured[i]==0) DAGmask[i]=1; else {DAGmask[i]=0;j++;}
	string DAGobsname = DAG->readObservationName(i).name;

	if (DAGobsname.find("->")!=string::npos) DAGobsname = DAGobsname.substr(DAGobsname.find("->")+2);

//if (allmeasured[i] && idnames[j] != DAGobsname)
	//		error("names do not match for observation at phenofile line (phe/geno) %i/+1 (%s/%s)\n",
	//			i+1,idnames[i].c_str(),DAGobsname.c_str());
	// fix thanks to Vadym Pinchuk
	if (allmeasured[i] && idnames[j] != DAGobsname)
	    error("names do not match for observation at phenofile line(phe/geno) %i/+1 (%s/%s)\n",
		  i+1,idnames[j].c_str(),DAGobsname.c_str());



    }
    nids = j+1;
    //fprintf(stdout,"in INI: %i %i\n",nids,npeople);
    if (nids != nmeasured) error("nids != mneasured (%i != %i)\n",nids,nmeasured);

}
void gendata::re_gendata(char * fname, int insnps, int ingpreds, int npeople,
			 int nmeasured, unsigned short int * allmeasured, int skipd,
			 std::string * idnames)
{
    nids = nmeasured;
    nsnps = insnps;
    ngpreds = ingpreds;
    DAG = NULL;
    //	int nids_all = npeople;

    G.reinit(nids,(nsnps*ngpreds));

    FILE * infile;

    if ((infile=fopen(fname,"r"))==NULL) {
	fprintf(stderr,"gendata: cannot open file %s\n",fname);
    }

    char tmp[100],tmpn[100];
    std::string tmpid,tmpstr;

    int k = 0;
    for (int i = 0;i<npeople;i++)
	if (allmeasured[i]==1)
	{
	    if (skipd>0)
	    {
		//				int ttt;
		char ttt[100];
		fscanf(infile,"%s", tmp);
		//				sscanf(tmp,"%d->%s",&ttt, tmpn);
		//		these changes are thanks to BMM & BP :)
		//				sscanf(tmp,"%s->%s",&ttt, tmpn);
		//				sscanf(tmp,"%[^->]->%[^->]",&ttt, tmpn);
		tmpstr = tmp;
		if (tmpstr.find("->")!=string::npos) {
		    sscanf(tmp,"%[^->]->%s",ttt, tmpn);
		    tmpid = tmpn;
		} else {
		    tmpid = tmpstr;
		    //fprintf(stdout,"%s;%s;%s;%s;%s\n",tmp,ttt,tmpn,tmpid.c_str(),idnames[k].c_str());
		}
		if (tmpid != idnames[k])
		{
		    fprintf(stderr,"phenofile and dosefile did not match at line %d ",i+2);
		    cerr << "(" << tmpid << " != " << idnames[k] << ")\n";
		    fclose(infile);
		    exit(1);
		}
	    }
	    for (int j=1;j<skipd;j++) {
		fscanf(infile,"%s", tmp);
	    }
	    for (int j=0;j<(nsnps*ngpreds);j++)
	    {
		int a = fscanf(infile,"%s", tmp);
		if (!a || a==EOF)
		{
		    fprintf(stderr,"cannot read dose-file: check skipd and ngpreds parameters\n");
		    fclose(infile);
		    exit(1);
		}
		G.put(atof(tmp),k,j);
	    }
	    k++;
	}
	else
	{
	    for (int j=0;j<skipd;j++) fscanf(infile,"%s", tmp);
	    for (int j=0;j<(nsnps*ngpreds);j++) fscanf(infile,"%s", tmp);
	}
    fclose(infile);
}
// HERE NEED A NEW CONSTRUCTOR BASED ON DATABELBASECPP OBJECT
gendata::~gendata()
{

    if (DAG != NULL) {delete DAG;delete [] DAGmask;}

    //		delete G;
}









class regdata
{
public:
    int nids;
    int ncov;
    int ngpreds;
    int noutcomes;
    unsigned short int * masked_data;
    mematrix<double> X;
    mematrix<double> Y;

    regdata()
    {
    }
    regdata(const regdata &obj)
    {
	nids = obj.nids;
	ncov = obj.ncov;
	ngpreds = obj.ngpreds;
	noutcomes = obj.noutcomes;
	X = obj.X;
	Y = obj.Y;
	masked_data = new unsigned short int [nids];
	for (int i=0;i<nids;i++) masked_data[i] = 0;
    }
    regdata(phedata &phed, gendata &gend, int snpnum)
    {
	nids = gend.nids;
	masked_data = new unsigned short int [nids];
	for (int i=0;i<nids;i++) masked_data[i] = 0;
	ngpreds = gend.ngpreds;
	if (snpnum>=0)
	    ncov = phed.ncov + ngpreds;
	else
	    ncov = phed.ncov;
	noutcomes = phed.noutcomes;
	X.reinit(nids,(ncov+1));
	Y.reinit(nids,noutcomes);
	for (int i=0;i<nids;i++)
	{
	    X.put(1.,i,0);
	    Y.put((phed.Y).get(i,0),i,0);
	}
	for (int j=1;j<=phed.ncov;j++)
	    for (int i=0;i<nids;i++)
		X.put((phed.X).get(i,j-1),i,j);
	if (snpnum>0)
	    for (int j=0;j<ngpreds;j++)
	    {
		float snpdata[nids];
		gend.get_var(snpnum*ngpreds+j,snpdata);
		for (int i=0;i<nids;i++)
		    X.put(snpdata[i],i,(ncov-ngpreds+1+j));
	    }
	//			for (int i=0;i<nids;i++)
	//				for (int j=0;j<ngpreds;j++)
	//					X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+1+j));
    }
    void update_snp(gendata &gend, int snpnum)
    {
	for (int j=0;j<ngpreds;j++)
	{
	    float snpdata[nids];
	    for (int i=0;i<nids;i++) masked_data[i]=0;
	    gend.get_var(snpnum*ngpreds+j,snpdata);
	    for (int i=0;i<nids;i++) {
		X.put(snpdata[i],i,(ncov+1-j-1));
		if (isnan(snpdata[i])) masked_data[i]=1;
	    }
	}
    }
    ~regdata()
    {
	delete [] masked_data;
	//		delete X;
	//		delete Y;
    }

    regdata get_unmasked_data()
    {
	regdata to; // = regdata(*this);
	int nmeasured = 0;
	for (int i=0;i<nids;i++)
	    if (masked_data[i]==0) nmeasured++;
	to.nids = nmeasured;
	//cout << to.nids << " in get_unmasked_data\n";
	to.ncov = ncov;
	to.ngpreds = ngpreds;
	to.noutcomes = noutcomes;
	int dim2Y = Y.ncol;
	int dim2X = X.ncol;
	(to.X).reinit(to.nids,dim2X);
	(to.Y).reinit(to.nids,dim2Y);

	int j = 0;
	for (int i=0;i<nids;i++)
	{
	    if (masked_data[i]==0) {
		for (int nc=0;nc<dim2X;nc++)
		    (to.X).put(X.get(i,nc),j,nc);
		for (int nc=0;nc<dim2Y;nc++)
		    (to.Y).put(Y.get(i,nc),j,nc);
		j++;
	    }
	}

	//delete [] to.masked_data;
	to.masked_data = new unsigned short int [to.nids];
	for (int i=0;i<to.nids;i++) to.masked_data[i] = 0;
	//fprintf(stdout,"get_unmasked: %i %i %i\n",to.nids,dim2X,dim2Y);
	return(to);
    }

    mematrix<double> extract_genotypes(void)
    {
	mematrix<double> out;
	out.reinit(X.nrow,ngpreds);
	for (int i=0;i<X.nrow;i++)
	    for (int j=0;j<ngpreds;j++)
		out[i*ngpreds+j] = X.get(i,(ncov-ngpreds+1+j));
	return out;
    }
};

// compare for sort of times
int cmpfun(const void *a, const void *b)
{
    double el1 = *(double*)a;
    double el2 = *(double*)b;
    if (el1>el2) return 1;
    if (el1<el2) return -1;
    if (el1==el2) return 0;
}






class coxph_data
{
public:
    int nids;
    int ncov;
    int ngpreds;
    mematrix<double> weights;
    mematrix<double> stime;
    mematrix<int>    sstat;
    mematrix<double> offset;
    mematrix<int>    strata;
    mematrix<double> X;
    mematrix<int>    order;
    unsigned short int * masked_data;

    coxph_data(){}

    coxph_data(const coxph_data &obj)
    {
	nids = obj.nids;
	ncov = obj.ncov;
	ngpreds = obj.ngpreds;
	weights = obj.weights;
	stime = obj.stime;
	sstat = obj.sstat;
	offset = obj.offset;
	strata = obj.strata;
	X = obj.X;
	order = obj.order;
	masked_data = new unsigned short int [nids];
	for (int i=0;i<nids;i++) masked_data[i] = 0;
    }
    coxph_data(phedata &phed, gendata &gend, int snpnum)
    {
	nids = gend.nids;
	masked_data = new unsigned short int [nids];
	for (int i=0;i<nids;i++) masked_data[i] = 0;
	ngpreds = gend.ngpreds;
	if (snpnum>=0)
	    ncov = phed.ncov + ngpreds;
	else
	    ncov = phed.ncov;
	if (phed.noutcomes != 2)
	{
	    fprintf(stderr,"coxph_data: number of outcomes should be 2 (now: %d)\n",phed.noutcomes);
	    exit(1);
	}
	//		X.reinit(nids,(ncov+1));
	X.reinit(nids,ncov);
	stime.reinit(nids,1);
	sstat.reinit(nids,1);
	weights.reinit(nids,1);
	offset.reinit(nids,1);
	strata.reinit(nids,1);
	order.reinit(nids,1);
	for (int i=0;i<nids;i++)
	{
	    //			X.put(1.,i,0);
	    stime[i] = (phed.Y).get(i,0);
	    sstat[i] = int((phed.Y).get(i,1));
	    if (sstat[i] != 1 & sstat[i]!=0)
	    {
		fprintf(stderr,"coxph_data: status not 0/1 (right order: id, fuptime, status ...)\n",phed.noutcomes);
		exit(1);
	    }
	}
	for (int j=0;j<phed.ncov;j++)
	    for (int i=0;i<nids;i++)
		X.put((phed.X).get(i,j),i,j);

	if (snpnum>0)
	    for (int j=0;j<ngpreds;j++)
	    {
		float snpdata[nids];
		gend.get_var(snpnum*ngpreds+j,snpdata);
		for (int i=0;i<nids;i++)
		    X.put(snpdata[i],i,(ncov-ngpreds+j));
	    }

	//			for (int i=0;i<nids;i++)
	//				for (int j=0;j<ngpreds;j++)
	//					X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+j));

	for (int i=0;i<nids;i++)
	{
	    weights[i] = 1.0;
	    offset[i] = 0.0;
	    strata[i] = 0;
	}
	// sort by time
	double tmptime[nids];
	int passed_sorted[nids];
	for (int i=0;i<nids;i++) {tmptime[i] = stime[i];passed_sorted[i]=0;}
	qsort(tmptime,nids,sizeof(double),cmpfun);
	for (int i=0;i<nids;i++)
	{
	    int passed = 0;
	    for (int j=0;j<nids;j++)
		if (tmptime[j] == stime[i])
		    if (!passed_sorted[j])
		    {
			order[i] = j;
			passed_sorted[j] = 1;
			passed = 1;
			break;
		    }
	    if (passed != 1)
	    {
		fprintf(stderr,"cannot recover element %d\n",i);
		exit(1);
	    }
	}
	stime = reorder(stime,order);
	sstat = reorder(sstat,order);
	weights = reorder(weights,order);
	strata = reorder(strata,order);
	offset = reorder(offset,order);
	X = reorder(X,order);
	X = transpose(X);
	//		X.print();
	//		offset.print();
	//		weights.print();
	//		stime.print();
	//		sstat.print();
    }
    void update_snp(gendata &gend, int snpnum)
    {
	// note this sorts by "order"!!!

	for (int j=0;j<ngpreds;j++)
	{
	    float snpdata[nids];
	    for (int i=0;i<nids;i++) masked_data[i]=0;
	    gend.get_var(snpnum*ngpreds+j,snpdata);
	    for (int i=0;i<nids;i++) {
		X.put(snpdata[i],(ncov-ngpreds+j),order[i]);
		if (isnan(snpdata[i])) masked_data[order[i]]=1;
	    }
	}
	//		for (int i=0;i<nids;i++)
	//			for (int j=0;j<ngpreds;j++)
	//				X.put((gend.G).get(i,(snpnum*ngpreds+j)),(ncov-ngpreds+j),order[i]);
    }
    ~coxph_data()
    {
	delete [] masked_data;
	//		delete X;
	//		delete sstat;
	//		delete stime;
	//		delete weights;
	//		delete offset;
	//		delete strata;
	//		delete order;
    }


    coxph_data get_unmasked_data()
    {
//		std::cout << " !!! in get_unmasked_data !!! ";
	coxph_data to; // = coxph_data(*this);
	// filter missing data

	int nmeasured = 0;
	for (int i=0;i<nids;i++)
	    if (masked_data[i]==0) nmeasured++;
	to.nids = nmeasured;
//		std::cout << "nmeasured=" << nmeasured << "\n";
	to.ncov = ncov;
//		std::cout << "ncov=" << ncov << "\n";
	to.ngpreds = ngpreds;
	int dim1X = X.nrow;
//		std::cout << "X.ncol=" << X.ncol << "\n";
//		std::cout << "X.nrow=" << X.nrow << "\n";
	(to.weights).reinit(to.nids,1);
	(to.stime).reinit(to.nids,1);
	(to.sstat).reinit(to.nids,1);
	(to.offset).reinit(to.nids,1);
	(to.strata).reinit(to.nids,1);
	(to.order).reinit(to.nids,1);
	(to.X).reinit(dim1X,to.nids);

//		std::cout << "(to.X).ncol=" << (to.X).ncol << "\n";
//		std::cout << "(to.X).nrow=" << (to.X).nrow << "\n";
//		std::cout << " !!! just before cycle !!! ";
	int j = 0;
	for (int i=0;i<nids;i++)
	{
//			std::cout << nids << " " << i << " " << masked_data[i] << "\n";
	    if (masked_data[i]==0) {
		(to.weights).put(weights.get(i,1),j,1);
		(to.stime).put(stime.get(i,1),j,1);
		(to.sstat).put(sstat.get(i,1),j,1);
		(to.offset).put(offset.get(i,1),j,1);
		(to.strata).put(strata.get(i,1),j,1);
		(to.order).put(order.get(i,1),j,1);
		for (int nc=0;nc<dim1X;nc++)
		    (to.X).put(X.get(nc,i),nc,j);
		j++;
	    }
	}
//		std::cout << " !!! just after cycle !!! ";

	//delete [] to.masked_data;
	to.masked_data = new unsigned short int [to.nids];
	for (int i=0;i<to.nids;i++) to.masked_data[i] = 0;
	//fprintf(stdout,"get_unmasked: %i %i %i\n",to.nids,dim2X,dim2Y);
	return(to);
    }
};






class mlinfo
{
public:
    int nsnps;
    std::string * name;
    std::string * A1;
    std::string * A2;
    double * Freq1;
    double * MAF;
    double * Quality;
    double * Rsq;
    std::string * map;

    mlinfo(char * filename, char * mapname)
    {
	FILE * infile = fopen(filename,"r");
	if (infile == NULL)
	{
	    fprintf(stderr,"mlinfo: cannot open file %s",filename);
	    exit(1);
	}
	char tmp[100];
	unsigned int nlin=0;
	while (fscanf(infile,"%s", tmp)!=EOF) {
	    nlin++;
	}
	fclose(infile);
	if (nlin % 7)
	{
	    fprintf(stderr,"mlinfo: number of columns != 7 in %s",filename);
	    exit(1);
	}
	nsnps = int(nlin/7) - 1;
	printf("Number of SNPs = %d\n",nsnps);
	name = new std::string [nsnps];
	A1 = new std::string [nsnps];
	A2 = new std::string [nsnps];
	Freq1 = new double [nsnps];
	MAF = new double [nsnps];
	Quality = new double [nsnps];
	Rsq = new double [nsnps];
	map = new std::string [nsnps];
	if ((infile = fopen(filename,"r"))==NULL)
	{
	    fprintf(stderr,"mlinfo: cannot open file %s",filename);
	    exit(1);
	}
	for (int i =0;i<7;i++) fscanf(infile,"%s", tmp);
	for (int i =0;i<nsnps;i++)
	{
	    fscanf(infile,"%s", tmp);
	    name[i] = tmp;
	    fscanf(infile,"%s", tmp);
	    A1[i] = tmp;
	    fscanf(infile,"%s", tmp);
	    A2[i] = tmp;
	    fscanf(infile,"%s", tmp);
	    Freq1[i] = atof(tmp);
	    fscanf(infile,"%s", tmp);
	    MAF[i] = atof(tmp);
	    fscanf(infile,"%s", tmp);
	    Quality[i] = atof(tmp);
	    fscanf(infile,"%s", tmp);
	    Rsq[i] = atof(tmp);
	    map[i] = "-999";
	}
	fclose(infile);
	if (mapname!=NULL)
	{
	    std::ifstream instr(mapname);
	    int BFS = 1000;
	    char line [BFS], tmp[BFS];
	    if (!instr.is_open())
	    {
		fprintf(stderr,"mlinfo: cannot open file %s",mapname);
		exit(1);
	    }
	    instr.getline(line,BFS);
	    for (int i=0;i<nsnps;i++)
	    {
		instr.getline(line,BFS);
		std::stringstream line_stream(line);
		line_stream >> tmp >> map[i];
	    }
	    instr.close();
	}
    }
    ~mlinfo()
    {
	delete [] name;
	delete [] A1;
	delete [] A2;
	delete [] Freq1;
	delete [] MAF;
	delete [] Quality;
	delete [] Rsq;
	delete [] map;
    }
};



//_________________________________________Maksim_start

class InvSigma
{

public:

    InvSigma(const char * filename_, phedata * phe)
    {
	filename = filename_;
	npeople  = phe->nids;
	std::ifstream myfile(filename_);
	char * line = new char[MAXIMUM_PEOPLE_AMOUNT];
	double val;
	std::string id;
	unsigned row=0, col=0;

	matrix.reinit(npeople,npeople);


	//idnames[k], if (allmeasured[i]==1)

	if (myfile.is_open())
	{
	    while(myfile.getline(line,MAXIMUM_PEOPLE_AMOUNT))
	    {



		std::stringstream line_stream(line);
		line_stream >> id;

		if(phe->idnames[row] != id) {std::cerr<<"error:in row "<<row<<" id="<<phe->idnames[row]<<" in inverse variance matrix but id="<<id<<" must be there. Wrong inverse variance matrix (only measured id must be there)\n";exit(1);}

		while (line_stream >> val)
		{
		    matrix.put(val, row, col);
		    col++;
		}

		if(col != npeople)
		{
		    fprintf(stderr,"error: inv file: Number of columns in row %d equals to %d but people amount is %d\n", row, col, npeople);
		    myfile.close();
		    exit(1);
		}
		col=0;
		row++;
	    }
	    myfile.close();
	}
	else
	{
	    fprintf(stderr,"error: inv file: cannot open file '%s'\n",filename_);
	}


    };


    ~InvSigma()
    {
	//af af
    }


    mematrix<double> & get_matrix(void)
    {
	return matrix;
    }


private:

    static const unsigned MAXIMUM_PEOPLE_AMOUNT = 1000000;
    unsigned npeople; //amount of people
    std::string filename;
    mematrix<double> matrix; //file is stored here

};
//________________________________________Maksim_end
