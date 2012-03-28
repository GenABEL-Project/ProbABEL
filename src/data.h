/*
 * data.h
 *
 *  Created on: Mar 8, 2012
 *      Author: mkooyman
 */

#ifndef DATA_H_
#define DATA_H_

extern bool is_interaction_excluded;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople);
#include "phedata.h"
#include "gendata.h"


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
	regdata(){}
	regdata(const regdata &obj) ;
	regdata(phedata &phed, gendata &gend, int snpnum);
	mematrix<double>  extract_genotypes();
	void update_snp(gendata &gend, int snpnum);
	regdata get_unmasked_data();
	~regdata();
private:

};

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
	coxph_data get_unmasked_data();
	coxph_data(){}
	coxph_data(const coxph_data &obj);
	coxph_data(phedata &phed, gendata &gend, int snpnum);
	void update_snp(gendata &gend, int snpnum);
	~coxph_data();

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
	mlinfo(){}
	mlinfo(char * filename, char * mapname);
	~mlinfo();
};

class InvSigma {

private:

	static const unsigned MAXIMUM_PEOPLE_AMOUNT = 1000000;
	unsigned npeople; //amount of people
	std::string filename;
	mematrix<double> matrix; //file is stored here
public:
	InvSigma(const char * filename_, phedata * phe);
	mematrix<double> & get_matrix();
	~InvSigma();
};

#endif /* DATA_H_ */
