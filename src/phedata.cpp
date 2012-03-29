#include <string>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cstdlib>
#include <phedata.h>

phedata::phedata(char * fname, int noutc, int npeople, int interaction,
        bool iscox) {
    setphedata(fname, noutc, npeople, interaction, iscox);
}

void phedata::set_is_interaction_excluded(bool int_exl) {
    is_interaction_excluded = int_exl;
}
void phedata::setphedata(char * fname, int noutc, int npeople, int interaction,
        bool iscox) {
    static const unsigned int BFS = 1000;
    std::ifstream myfile(fname);
    char line[BFS];
    char tmp[100];
    noutcomes = noutc;

    int nphenocols = 0;
    int savenpeople = npeople;
    npeople = 0;
    if (myfile.is_open()) {
        myfile.getline(line, BFS);
        std::stringstream line_stream(line);
        //			std::cout << line << "\n ";
        while (line_stream >> tmp) {

            nphenocols++;
            //				std::cout << tmp << " " << nphenocols << " ";
        }
        while (myfile.getline(line, BFS)) {
            int tmplins = 0;
            std::stringstream line_stream(line);
            while (line_stream >> tmp)
                tmplins++;
            if (tmplins != nphenocols) {
                std::cerr << "phenofile: number of variables different from "
                        << nphenocols << " in line " << tmplins << endl;
                myfile.close();
                exit(1);
            }
            npeople++;
        };
        myfile.close();
    } else {
        std::cerr << "Unable to open file " << fname << endl;
        exit(1);
    }
    std::cout << "Actual number of people in phenofile = " << npeople;
    if (savenpeople > 0) {
        npeople = savenpeople;
        std::cout << "; using only " << npeople << " first\n";
    } else {
        std::cout << "; using all of these\n";
    }

    ncov = nphenocols - 1 - noutcomes;
    nids_all = npeople;
    model_terms = new std::string[ncov + 2];

    // first pass -- find unmeasured people
    std::ifstream infile(fname);
    if (!infile) {
        std::cerr << "phedata: cannot open file " << fname << endl;
    }

    infile >> tmp;
    model = "( ";
    infile >> tmp;
    model = model + tmp;
    for (int i = 1; i < noutcomes; i++) {
        infile >> tmp;
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

    if (nphenocols > noutcomes + 1) {
        infile >> tmp;
        model = model + tmp;
        model_terms[n_model_terms++] = tmp;
        for (int i = (2 + noutcomes); i < nphenocols; i++) {
            infile >> tmp;

            //				if(iscox && ) {if(n_model_terms+1 == interaction-1) {continue;} }
            //				else      {if(n_model_terms+1 == interaction) {continue;} }
            model = model + " + ";
            model = model + tmp;
            model_terms[n_model_terms++] = tmp;
        }
    }
    model = model + " + SNP_A1";
    if (interaction != 0) {
        if (iscox) {
            model = model + " + " + model_terms[interaction - 1] + "*SNP_A1";
        } else {
            model = model + " + " + model_terms[interaction] + "*SNP_A1";
        }
    }
    model_terms[n_model_terms++] = "SNP_A1";

    if (is_interaction_excluded) // exclude covariates from covariate names
    {
        if (iscox) {
            std::cout << "model is running without "
                    << model_terms[interaction - 1] << ", term\n";
        } else {
            std::cout << "model is running without " << model_terms[interaction]
                    << ", term\n";
        }
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

    allmeasured = new unsigned short int[npeople];
    nids = 0;
    for (int i = 0; i < npeople; i++) {
        allmeasured[i] = 1;
        for (int j = 0; j < nphenocols; j++) {
            infile >> tmp;
            if (j > 0 && (tmp[0] == 'N' || tmp[0] == 'n'))
                allmeasured[i] = 0;
        }
        if (allmeasured[i] == 1)
            nids++;
    }
    infile.close();
    //		printf("npeople = %d, no. all measured = %d\n",nids_all,nids);

    // allocate objects
    int ntmpcov = 1;
    if (ncov > 0)
        ntmpcov = ncov;
    idnames = new std::string[nids];
    X.reinit(nids, ntmpcov);
    Y.reinit(nids, noutcomes);

    // second pass -- read the data
    infile.open(fname);
    if (!infile) {
        std::cerr << "phedata: cannot open file " << fname << endl;
        exit(1);
    }

    for (int i = 0; i < nphenocols; i++) {
        infile >> tmp;
    }

    int k = 0;
    int m = 0;
    for (int i = 0; i < npeople; i++)
        if (allmeasured[i] == 1) {
            infile >> tmp;
            idnames[m] = tmp;
            for (int j = 0; j < noutcomes; j++) {
                infile >> tmp;
                Y.put(atof(tmp), m, j);
            }
            for (int j = (1 + noutcomes); j < nphenocols; j++) {
                infile >> tmp;
                X.put(atof(tmp), m, (j - 1 - noutcomes));
            }
            m++;
        } else
            for (int j = 0; j < nphenocols; j++)
                infile >> tmp;
    infile.close();
}
phedata::~phedata() {
    //		delete X;
    //		delete Y;
    //		delete [] allmeasured;
}

