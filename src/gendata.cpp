/*
 * gendata.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: mkooyman
 */
#include <string>
#include "gendata.h"
#include "fvlib/FileVector.h"
#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif
#include "utilities.h"

void gendata::get_var(int var, float * data)
{
    if (DAG == NULL)
        for (int i = 0; i < G.nrow; i++)
            data[i] = G.get(i, var);
    else if (DAG != NULL)
    {
        float tmpdata[DAG->getNumObservations()];
        DAG->readVariableAs((unsigned long int) var, tmpdata);
        unsigned int j = 0;
        for (unsigned int i = 0; i < DAG->getNumObservations(); i++)
            if (!DAGmask[i])
                data[j++] = tmpdata[i];
        // std::cout << j << " " << DAG->get_nobservations() << " "
        //           << nids << "\n";
    }
    else
    {
        report_error("cannot get gendata");
    }
}

gendata::gendata() : nsnps(0), nids(0), ngpreds(0), DAG(NULL), DAGmask(NULL)
{
}

void gendata::re_gendata(string filename, unsigned int insnps,
                         unsigned int ingpreds, unsigned int npeople,
                         unsigned int nmeasured,
                         unsigned short int * allmeasured,
                         std::string * idnames)
{
    nsnps = insnps;
    ngpreds = ingpreds;
    DAG = new FileVector(filename, 128, true);
    DAGmask = new unsigned short int[DAG->getNumObservations()];
    if (DAG->getNumObservations() != npeople)
        report_error("dimension of fvf-data and phenotype data do not match\n");

    if (DAG->getNumVariables() != insnps * ingpreds)
        report_error("dimension of fvf-data and mlinfo data do not match\n");

    long int j = -1;

    for (unsigned int i = 0; i < npeople; i++)
    {
        if (allmeasured[i] == 0)
            DAGmask[i] = 1;
        else
        {
            DAGmask[i] = 0;
            j++;
        }
        string DAGobsname = DAG->readObservationName(i).name;

        if (DAGobsname.find("->") != string::npos)
            DAGobsname = DAGobsname.substr(DAGobsname.find("->") + 2);

        // if (allmeasured[i] && idnames[j] != DAGobsname)
        //  std::cerr << "names do not match for observation at phenofile "
        //            << "line (phe/geno) " << i+1 << "/+1 ("
        //            << idnames[i].c_str() << "/"
        //            << DAGobsname.c_str() << ")\n";
        // fix thanks to Vadym Pinchuk
        if (allmeasured[i] && idnames[j] != DAGobsname)
            report_error(
                "names do not match for observation at phenofile line(phe/geno) %i/+1 (%s/%s)\n",
                i + 1, idnames[j].c_str(), DAGobsname.c_str());

    }
    nids = j + 1;
    // std::cout << "in INI: " << nids << " " << npeople << "\n";
    if (nids != nmeasured)
        report_error("nids != mneasured (%i != %i)\n", nids, nmeasured);
}

void gendata::re_gendata(char * fname, unsigned int insnps,
                         unsigned int ingpreds, unsigned int npeople,
                         unsigned int nmeasured,
                         unsigned short int * allmeasured, int skipd,
                         std::string * idnames)
{
    nids    = nmeasured;
    nsnps   = insnps;
    ngpreds = ingpreds;
    DAG     = NULL;
    //	int nids_all = npeople;

    G.reinit(nids, (nsnps * ngpreds));

    std::ifstream infile;

    infile.open(fname);
    if (!infile)
    {
        std::cerr << "gendata: cannot open file " << fname << endl;
    }

    std::string tmpid, tmpstr;

    int k = 0;
    for (unsigned int i = 0; i < npeople; i++)
    {
        if (allmeasured[i] == 1)
        {
            if (skipd > 0)
            {
                // Read the genotype data and look for the signature
                // arrow of MaCH/minimac. If found only use the part
                // after the arrow as ID.
                infile >> tmpstr;
                size_t strpos = tmpstr.find("->");
                if (strpos != string::npos)
                {
                    tmpid = tmpstr.substr(strpos+2, string::npos);
                }
                else
                {
                    tmpid = tmpstr;
                }
                if (tmpid != idnames[k])
                {
                    cerr << "phenotype file and dose or probability file "
                         << "did not match at line " << i + 2 << "(" << tmpid
                         << " != " << idnames[k] << ")" << endl;
                    infile.close();
                    exit(1);
                }
            }

            for (int j = 1; j < skipd; j++)
            {
                infile >> tmpstr;
            }

            for (int j = 0; j < (nsnps * ngpreds); j++)
            {
                if (infile.good())
                {
                    infile >> tmpstr;
                    // tmpstr contains the dosage in string form. Convert
                    // it to float (if tmpstr is NaN it will be set to nan).
                    // Note that Valgrind 3.7.0 gives "Invalid read of
                    // size 8" error messages here. A bug in Valgrind?!
                    float dosage = strtod(tmpstr.c_str(), (char **) NULL);
                    G.put(dosage, k, j);
                }
                else
                {
                    std::cerr << "cannot read dose-file: "
                              << "check skipd and ngpreds parameters\n";
                    infile.close();
                    exit(1);
                }
            }
            k++;
        }
        else
        {
            for (int j = 0; j < skipd; j++)
                infile >> tmpstr;
            for (int j = 0; j < (nsnps * ngpreds); j++)
                infile >> tmpstr;
        }
    }
    infile.close();
}
// HERE NEED A NEW CONSTRUCTOR BASED ON DATABELBASECPP OBJECT
gendata::~gendata()
{
    if (DAG != NULL)
    {
        delete DAG;
        delete[] DAGmask;
    }

    //		delete G;
}