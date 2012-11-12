#include <string>
#include <sstream>
#include <fstream>

#include "fvlib/AbstractMatrix.h"
#include "fvlib/CastUtils.h"
#include "fvlib/const.h"
#include "fvlib/convert_util.h"
#include "fvlib/FileVector.h"
#include "fvlib/frutil.h"
#include "fvlib/frversion.h"
#include "fvlib/Logger.h"
#include "fvlib/Transposer.h"
#include "phedata.h"
#include "gendata.h"
#include "data.h"

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif
#include "utilities.h"
using namespace std;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople)
{
//TODO: unused variables remove them for good if there is no reason to keep them
//int ncov = nphenocols - 2;
//int nids_all = npeople;

// first pass -- find unmeasured people
    std::ifstream infile(fname);
    if (!infile)
    {
        std::cerr << "Nmeasured: cannot open file " << fname << endl;
    }
    char tmp[100];

    for (int i = 0; i < nphenocols; i++)
    {
        infile >> tmp;
    }

    unsigned short int * allmeasured = new unsigned short int[npeople];
    int nids = 0;
    for (int i = 0; i < npeople; i++)
    {
        allmeasured[i] = 1;
        infile >> tmp;
        for (int j = 1; j < nphenocols; j++)
        {
            infile >> tmp;
            if (tmp[0] == 'N' || tmp[0] == 'n')
                allmeasured[i] = 0;
        }
        if (allmeasured[i] == 1)
            nids++;
    }
    infile.close();

    delete[] allmeasured;

    return (nids);
}

mlinfo::mlinfo(char * filename, char * mapname)
{

    char tmp[100];
    unsigned int nlin = 0;
    std::ifstream infile(filename);
    if (infile.is_open())
    {
        while (infile.good())
        {
            infile >> tmp;
            nlin++;
        }
        nlin--; // Subtract one, the previous loop added 1 too much
    } else
    {
        std::cerr << "mlinfo: cannot open file " << filename << endl;
        exit(1);
    }
    infile.close();

    if (nlin % 7)
    {
        std::cerr << "mlinfo: number of columns != 7 in " << filename << endl;
        exit(1);
    }
    nsnps = int(nlin / 7) - 1;
    std::cout << "Number of SNPs = " << nsnps << endl;
    name = new std::string[nsnps];
    A1 = new std::string[nsnps];
    A2 = new std::string[nsnps];
    Freq1 = new double[nsnps];
    MAF = new double[nsnps];
    Quality = new double[nsnps];
    Rsq = new double[nsnps];
    map = new std::string[nsnps];

    infile.open(filename);
    if (!infile)
    { // file couldn't be opened
        std::cerr << "mlinfo: cannot open file " << filename << endl;
        exit(1);
    }
    /* Read the header and discard it */
    for (int i = 0; i < 7; i++)
        infile >> tmp;

    for (int i = 0; i < nsnps; i++)
    {
        infile >> tmp;
        name[i] = tmp;
        infile >> tmp;
        A1[i] = tmp;
        infile >> tmp;
        A2[i] = tmp;
        infile >> tmp;
        Freq1[i] = atof(tmp);
        infile >> tmp;
        MAF[i] = atof(tmp);
        infile >> tmp;
        Quality[i] = atof(tmp);
        infile >> tmp;
        Rsq[i] = atof(tmp);
        map[i] = "-999";
    }
    infile.close();

    if (mapname != NULL)
    {
        std::ifstream instr(mapname);
        int BFS = 1000;
        char line[BFS], tmp[BFS];
        if (!instr.is_open())
        {
            std::cerr << "mlinfo: cannot open file " << mapname << endl;
            exit(1);
        }
        instr.getline(line, BFS);
        for (int i = 0; i < nsnps; i++)
        {
            instr.getline(line, BFS);
            std::stringstream line_stream(line);
            line_stream >> tmp >> map[i];
        }
        instr.close();
    }
}
mlinfo::~mlinfo()
{
    delete[] mlinfo::name;
    delete[] mlinfo::A1;
    delete[] mlinfo::A2;
    delete[] mlinfo::Freq1;
    delete[] mlinfo::MAF;
    delete[] mlinfo::Quality;
    delete[] mlinfo::Rsq;
    delete[] mlinfo::map;
}

//_________________________________________Maksim_start

InvSigma::InvSigma(const char * filename_, phedata * phe)
{
    filename = filename_;
    npeople = phe->nids;
    std::ifstream myfile(filename_);
    char * line = new char[MAXIMUM_PEOPLE_AMOUNT];
    std::string id;

    matrix.reinit(npeople, npeople);

//idnames[k], if (allmeasured[i]==1)

    if (myfile.is_open())
    {
        double val;
        unsigned row = 0;
        while (myfile.getline(line, MAXIMUM_PEOPLE_AMOUNT))
        {

            std::stringstream line_stream(line);
            line_stream >> id;

            if (phe->idnames[row] != id)
            {
                std::cerr << "error:in row " << row << " id="
                        << phe->idnames[row]
                        << " in inverse variance matrix but id=" << id
                        << " must be there. Wrong inverse variance matrix (only measured id must be there)\n";
                exit(1);
            }
            unsigned col = 0;
            while (line_stream >> val)
            {
                matrix.put(val, row, col);
                col++;
            }

            if (col != npeople)
            {
                fprintf(stderr,
                        "error: inv file: Number of columns in row %d equals to %d but people amount is %d\n",
                        row, col, npeople);
                myfile.close();
                exit(1);
            }
            col = 0;
            row++;
        }
        myfile.close();
    } else
    {
        fprintf(stderr, "error: inv file: cannot open file '%s'\n", filename_);
    }

    delete[] line;
}
;

InvSigma::~InvSigma()
{
//af af
}

mematrix<double> & InvSigma::get_matrix(void)
{
    return matrix;
}

//________________________________________Maksim_end
