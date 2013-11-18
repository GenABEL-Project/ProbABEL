/**
 * \file   data.h
 * \author mkooyman
 *
 * \brief Contains several classes we didn't put somewhere else yet
 */

#ifndef DATA_H_
#define DATA_H_
#include <string>

extern bool is_interaction_excluded;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople);
#include "phedata.h"
#include "gendata.h"

/**
 * \brief Data from the mlinfo file.
 *
 */
class mlinfo {
public:
    int nsnps;                  /**< Number of SNPs */
    std::string * name;         /**< Array of SNP names */
    std::string * A1;           /**< Array with the first allele */
    std::string * A2;           /**< Array with the second allele */
    double * Freq1;
    double * MAF;               /**< The minor allele frequency */
    double * Quality;           /**< The imputation quality metric */
    double * Rsq;               /**< The imputation \f$R^2\f$ */
    std::string * map;          /**< Array with the SNP positions */
    mlinfo()
    {
        Freq1 = NULL;
        MAF = NULL;
        Quality = NULL;
        Rsq = NULL;
        nsnps = 0;
        A1 = NULL;
        A2 = NULL;
        name = NULL;
        map = NULL;
    }
    mlinfo(char * filename, char * mapname);
    ~mlinfo();
};

class InvSigma {
private:
    static const unsigned MAXIMUM_PEOPLE_AMOUNT = 1000000;
    unsigned int npeople;       /* number of people */
    std::string filename;
    mematrix<double> matrix;    /* file is stored here */

public:
    InvSigma(const char * filename_, phedata * phe);
    mematrix<double> & get_matrix();
    ~InvSigma();
};

#endif /* DATA_H_ */
