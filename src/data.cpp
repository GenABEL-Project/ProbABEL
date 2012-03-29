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

#include "mematrix.h"
#include "mematri1.h"
#include "utilities.h"
using namespace std;

extern bool is_interaction_excluded;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople) {
    int ncov = nphenocols - 2;
    int nids_all = npeople;

    // first pass -- find unmeasured people
    std::ifstream infile(fname);
    if (!infile) {
        std::cerr << "Nmeasured: cannot open file " << fname << endl;
    }
    char tmp[100];

    for (int i = 0; i < nphenocols; i++) {
        infile >> tmp;
    }

    unsigned short int * allmeasured = new unsigned short int[npeople];
    int nids = 0;
    for (int i = 0; i < npeople; i++) {
        allmeasured[i] = 1;
        infile >> tmp;
        for (int j = 1; j < nphenocols; j++) {
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

regdata::regdata(const regdata &obj) {
    nids = obj.nids;
    ncov = obj.ncov;
    ngpreds = obj.ngpreds;
    noutcomes = obj.noutcomes;
    X = obj.X;
    Y = obj.Y;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++) {
        masked_data[i] = 0;
    }
}
regdata::regdata(phedata &phed, gendata &gend, int snpnum) {
    nids = gend.nids;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++) {
        masked_data[i] = 0;
    }

    ngpreds = gend.ngpreds;
    if (snpnum >= 0) {
        ncov = phed.ncov + ngpreds;
    } else {
        ncov = phed.ncov;
    }
    noutcomes = phed.noutcomes;
    X.reinit(nids, (ncov + 1));
    Y.reinit(nids, noutcomes);
    for (int i = 0; i < nids; i++) {
        X.put(1., i, 0);
        Y.put((phed.Y).get(i, 0), i, 0);
    }
    for (int j = 1; j <= phed.ncov; j++)
        for (int i = 0; i < nids; i++)
            X.put((phed.X).get(i, j - 1), i, j);
    if (snpnum > 0)
        for (int j = 0; j < ngpreds; j++) {
            float snpdata[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
                X.put(snpdata[i], i, (ncov - ngpreds + 1 + j));
        }
    //			for (int i=0;i<nids;i++)
    //				for (int j=0;j<ngpreds;j++)
    //					X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+1+j));
}
void regdata::update_snp(gendata &gend, int snpnum) {
    for (int j = 0; j < ngpreds; j++) {
        float snpdata[nids];
        for (int i = 0; i < nids; i++)
            masked_data[i] = 0;
        gend.get_var(snpnum * ngpreds + j, snpdata);
        for (int i = 0; i < nids; i++) {
            X.put(snpdata[i], i, (ncov + 1 - j - 1));
            if (isnan(snpdata[i]))
                masked_data[i] = 1;
        }
    }
}
regdata::~regdata() {
    delete[] regdata::masked_data;
    //		delete X;
    //		delete Y;
}

regdata regdata::get_unmasked_data() {
    regdata to; // = regdata(*this);
    int nmeasured = 0;
    for (int i = 0; i < nids; i++)
        if (masked_data[i] == 0)
            nmeasured++;
    to.nids = nmeasured;
    //cout << to.nids << " in get_unmasked_data\n";
    to.ncov = ncov;
    to.ngpreds = ngpreds;
    to.noutcomes = noutcomes;
    int dim2Y = Y.ncol;
    int dim2X = X.ncol;
    (to.X).reinit(to.nids, dim2X);
    (to.Y).reinit(to.nids, dim2Y);

    int j = 0;
    for (int i = 0; i < nids; i++) {
        if (masked_data[i] == 0) {
            for (int nc = 0; nc < dim2X; nc++)
                (to.X).put(X.get(i, nc), j, nc);
            for (int nc = 0; nc < dim2Y; nc++)
                (to.Y).put(Y.get(i, nc), j, nc);
            j++;
        }
    }

    //delete [] to.masked_data;
    to.masked_data = new unsigned short int[to.nids];
    for (int i = 0; i < to.nids; i++)
        to.masked_data[i] = 0;
    //fprintf(stdout,"get_unmasked: %i %i %i\n",to.nids,dim2X,dim2Y);
    return (to);
}

mematrix<double> regdata::extract_genotypes(void) {
    mematrix<double> out;
    out.reinit(X.nrow, ngpreds);
    for (int i = 0; i < X.nrow; i++)
        for (int j = 0; j < ngpreds; j++)
            out[i * ngpreds + j] = X.get(i, (ncov - ngpreds + 1 + j));
    return out;
}

// compare for sort of times
int cmpfun(const void *a, const void *b) {
    double el1 = *(double*) a;
    double el2 = *(double*) b;
    if (el1 > el2)
        return 1;
    if (el1 < el2)
        return -1;
    if (el1 == el2)
        return 0;
}

coxph_data::coxph_data(const coxph_data &obj) {
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
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++)
        masked_data[i] = 0;
}
coxph_data::coxph_data(phedata &phed, gendata &gend, int snpnum) {
    nids = gend.nids;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++)
        masked_data[i] = 0;
    ngpreds = gend.ngpreds;
    if (snpnum >= 0)
        ncov = phed.ncov + ngpreds;
    else
        ncov = phed.ncov;
    if (phed.noutcomes != 2) {
        fprintf(stderr,
                "coxph_data: number of outcomes should be 2 (now: %d)\n",
                phed.noutcomes);
        exit(1);
    }
    //		X.reinit(nids,(ncov+1));
    X.reinit(nids, ncov);
    stime.reinit(nids, 1);
    sstat.reinit(nids, 1);
    weights.reinit(nids, 1);
    offset.reinit(nids, 1);
    strata.reinit(nids, 1);
    order.reinit(nids, 1);
    for (int i = 0; i < nids; i++) {
        //			X.put(1.,i,0);
        stime[i] = (phed.Y).get(i, 0);
        sstat[i] = int((phed.Y).get(i, 1));
        if (sstat[i] != 1 & sstat[i] != 0) {
            fprintf(stderr,
                    "coxph_data: status not 0/1 (right order: id, fuptime, status ...)\n",
                    phed.noutcomes);
            exit(1);
        }
    }
    for (int j = 0; j < phed.ncov; j++)
        for (int i = 0; i < nids; i++)
            X.put((phed.X).get(i, j), i, j);

    if (snpnum > 0)
        for (int j = 0; j < ngpreds; j++) {
            float snpdata[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
                X.put(snpdata[i], i, (ncov - ngpreds + j));
        }

    //			for (int i=0;i<nids;i++)
    //				for (int j=0;j<ngpreds;j++)
    //					X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+j));

    for (int i = 0; i < nids; i++) {
        weights[i] = 1.0;
        offset[i] = 0.0;
        strata[i] = 0;
    }
    // sort by time
    double tmptime[nids];
    int passed_sorted[nids];
    for (int i = 0; i < nids; i++) {
        tmptime[i] = stime[i];
        passed_sorted[i] = 0;
    }
    qsort(tmptime, nids, sizeof(double), cmpfun);
    for (int i = 0; i < nids; i++) {
        int passed = 0;
        for (int j = 0; j < nids; j++)
            if (tmptime[j] == stime[i])
                if (!passed_sorted[j]) {
                    order[i] = j;
                    passed_sorted[j] = 1;
                    passed = 1;
                    break;
                }
        if (passed != 1) {
            fprintf(stderr, "cannot recover element %d\n", i);
            exit(1);
        }
    }
    stime = reorder(stime, order);
    sstat = reorder(sstat, order);
    weights = reorder(weights, order);
    strata = reorder(strata, order);
    offset = reorder(offset, order);
    X = reorder(X, order);
    X = transpose(X);
    //		X.print();
    //		offset.print();
    //		weights.print();
    //		stime.print();
    //		sstat.print();
}
void coxph_data::update_snp(gendata &gend, int snpnum) {
    // note this sorts by "order"!!!

    for (int j = 0; j < ngpreds; j++) {
        float snpdata[nids];
        for (int i = 0; i < nids; i++)
            masked_data[i] = 0;
        gend.get_var(snpnum * ngpreds + j, snpdata);
        for (int i = 0; i < nids; i++) {
            X.put(snpdata[i], (ncov - ngpreds + j), order[i]);
            if (isnan(snpdata[i]))
                masked_data[order[i]] = 1;
        }
    }
    //		for (int i=0;i<nids;i++)
    //			for (int j=0;j<ngpreds;j++)
    //				X.put((gend.G).get(i,(snpnum*ngpreds+j)),(ncov-ngpreds+j),order[i]);
}
coxph_data::~coxph_data() {
    delete[] coxph_data::masked_data;
    //		delete X;
    //		delete sstat;
    //		delete stime;
    //		delete weights;
    //		delete offset;
    //		delete strata;
    //		delete order;
}

coxph_data coxph_data::get_unmasked_data() {
//		std::cout << " !!! in get_unmasked_data !!! ";
    coxph_data to; // = coxph_data(*this);
    // filter missing data

    int nmeasured = 0;
    for (int i = 0; i < nids; i++)
        if (masked_data[i] == 0)
            nmeasured++;
    to.nids = nmeasured;
//		std::cout << "nmeasured=" << nmeasured << "\n";
    to.ncov = ncov;
//		std::cout << "ncov=" << ncov << "\n";
    to.ngpreds = ngpreds;
    int dim1X = X.nrow;
//		std::cout << "X.ncol=" << X.ncol << "\n";
//		std::cout << "X.nrow=" << X.nrow << "\n";
    (to.weights).reinit(to.nids, 1);
    (to.stime).reinit(to.nids, 1);
    (to.sstat).reinit(to.nids, 1);
    (to.offset).reinit(to.nids, 1);
    (to.strata).reinit(to.nids, 1);
    (to.order).reinit(to.nids, 1);
    (to.X).reinit(dim1X, to.nids);

//		std::cout << "(to.X).ncol=" << (to.X).ncol << "\n";
//		std::cout << "(to.X).nrow=" << (to.X).nrow << "\n";
//		std::cout << " !!! just before cycle !!! ";
    int j = 0;
    for (int i = 0; i < nids; i++) {
//			std::cout << nids << " " << i << " " << masked_data[i] << "\n";
        if (masked_data[i] == 0) {
            (to.weights).put(weights.get(i, 1), j, 1);
            (to.stime).put(stime.get(i, 1), j, 1);
            (to.sstat).put(sstat.get(i, 1), j, 1);
            (to.offset).put(offset.get(i, 1), j, 1);
            (to.strata).put(strata.get(i, 1), j, 1);
            (to.order).put(order.get(i, 1), j, 1);
            for (int nc = 0; nc < dim1X; nc++)
                (to.X).put(X.get(nc, i), nc, j);
            j++;
        }
    }
//		std::cout << " !!! just after cycle !!! ";

    //delete [] to.masked_data;
    to.masked_data = new unsigned short int[to.nids];
    for (int i = 0; i < to.nids; i++)
        to.masked_data[i] = 0;
    //fprintf(stdout,"get_unmasked: %i %i %i\n",to.nids,dim2X,dim2Y);
    return (to);
}

mlinfo::mlinfo(char * filename, char * mapname) {

    char tmp[100];
    unsigned int nlin = 0;
    std::ifstream infile(filename);
    if (infile.is_open()) {
        while (infile.good()) {
            infile >> tmp;
            nlin++;
        }
        nlin--; // Subtract one, the previous loop added 1 too much
    } else {
        std::cerr << "mlinfo: cannot open file " << filename << endl;
        exit(1);
    }
    infile.close();

    if (nlin % 7) {
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
    if (!infile) { // file couldn't be opened
        std::cerr << "mlinfo: cannot open file " << filename << endl;
        exit(1);
    }
    /* Read the header and discard it */
    for (int i = 0; i < 7; i++)
        infile >> tmp;

    for (int i = 0; i < nsnps; i++) {
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

    if (mapname != NULL) {
        std::ifstream instr(mapname);
        int BFS = 1000;
        char line[BFS], tmp[BFS];
        if (!instr.is_open()) {
            std::cerr << "mlinfo: cannot open file " << mapname << endl;
            exit(1);
        }
        instr.getline(line, BFS);
        for (int i = 0; i < nsnps; i++) {
            instr.getline(line, BFS);
            std::stringstream line_stream(line);
            line_stream >> tmp >> map[i];
        }
        instr.close();
    }
}
mlinfo::~mlinfo() {
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

InvSigma::InvSigma(const char * filename_, phedata * phe) {
filename = filename_;
npeople = phe->nids;
std::ifstream myfile(filename_);
char * line = new char[MAXIMUM_PEOPLE_AMOUNT];
double val;
std::string id;
unsigned row = 0, col = 0;

matrix.reinit(npeople, npeople);

//idnames[k], if (allmeasured[i]==1)

if (myfile.is_open()) {
    while (myfile.getline(line, MAXIMUM_PEOPLE_AMOUNT)) {

        std::stringstream line_stream(line);
        line_stream >> id;

        if (phe->idnames[row] != id) {
            std::cerr << "error:in row " << row << " id=" << phe->idnames[row]
                    << " in inverse variance matrix but id=" << id
                    << " must be there. Wrong inverse variance matrix (only measured id must be there)\n";
            exit(1);
        }

        while (line_stream >> val) {
            matrix.put(val, row, col);
            col++;
        }

        if (col != npeople) {
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
} else {
    fprintf(stderr, "error: inv file: cannot open file '%s'\n", filename_);
}

delete[] line;
}
;

InvSigma::~InvSigma() {
//af af
}

mematrix<double> & InvSigma::get_matrix(void) {
return matrix;
}

//________________________________________Maksim_end
