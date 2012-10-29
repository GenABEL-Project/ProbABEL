/*
 * comand_line_settings.cpp
 *
 *  Created on: Apr 1, 2012
 *      Author: mkooyman
 */

#include <getopt.h>
#include <string>
#include <iostream>
#include "usage.h"
#include "comand_line_settings.h"

// config.h and fvlib/FileVector.h are included for the upper case variables
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "fvlib/FileVector.h"

using namespace std;

string cmdvars::getStrGenfilename() const
{
    return str_genfilename;
}

int cmdvars::getAllcov() const
{
    return allcov;
}

string cmdvars::getChrom() const
{
    return chrom;
}

char* cmdvars::getGenfilename() const
{
    return genfilename;
}

int cmdvars::getInteraction() const
{
    return interaction;
}

char* cmdvars::getInverseFilename() const
{
    return inverse_filename;
}

bool cmdvars::isIscox() const
{
    return iscox;
}

int cmdvars::getIsFvf() const
{
    return isFVF;
}

char* cmdvars::getMapfilename() const
{
    return mapfilename;
}

char* cmdvars::getMlinfofilename() const
{
    return mlinfofilename;
}


int cmdvars::getNgpreds() const
{
    return ngpreds;
}

int cmdvars::getNohead() const
{
    return nohead;
}

int cmdvars::getNoutcomes() const
{
    return noutcomes;
}

int cmdvars::getNpeople() const
{
    return npeople;
}

char* cmdvars::getOutfilename() const
{
    return outfilename;
}

char* cmdvars::getProgramName() const
{
    return program_name;
}

int cmdvars::getRobust() const
{
    return robust;
}

int cmdvars::getScore() const
{
    return score;
}

string cmdvars::getSep() const
{
    return sep;
}

int cmdvars::getSkipd() const
{
    return skipd;
}

char* cmdvars::getPhefilename() const
{
    return phefilename;
}

void cmdvars::set_variables(int argc, char * argv[])
{
    int next_option;
    const char * const short_options = "p:i:d:m:n:c:o:s:t:g:a:erlh:b:vu";
    //b - interaction parameter
    // ADD --fv FLAG (FILEVECTOR), IN WHICH CASE USE ALTERNATIVE
    // CONSTRUCTOR FOR GENDATA
    const struct option long_options[] =
    {
    { "pheno", 1, NULL, 'p' },
    { "info", 1, NULL, 'i' },
    { "dose", 1, NULL, 'd' },
    { "map", 1, NULL, 'm' },
    { "nids", 1, NULL, 'n' },
    { "chrom", 1, NULL, 'c' },
    { "out", 1, NULL, 'o' },
    { "skipd", 1, NULL, 's' },
    { "ntraits", 1, NULL, 't' },
    { "ngpreds", 1, NULL, 'g' },
    { "separat", 1, NULL, 'a' },
    { "score", 0, NULL, 'r' },
    { "no-head", 0, NULL, 'e' },
    { "allcov", 0, NULL, 'l' },
    { "help", 0, NULL, 'h' },
    { "interaction", 1, NULL, 'b' },
    { "interaction_only", 1, NULL, 'k' },
    { "mmscore", 1, NULL, 'v' },
    { "robust", 0, NULL, 'u' },
    { NULL, 0, NULL, 0 } };
    char * program_name = argv[0];
    fprintf(stdout, "Usage: %s options\n", PACKAGE_VERSION);

    do
    {
        next_option = getopt_long(argc, argv, short_options, long_options,
                NULL);

        switch (next_option)
        {
        case 'h':
            print_help(program_name, 0);
        case 'p':
            phefilename = optarg;
            neco[0] = 1;
            fprintf(stdout, "phenoint\n");
            break;
        case 'i':
            mlinfofilename = optarg;
            neco[1] = 1;
            break;
        case 'd':
            genfilename = optarg;
            neco[2] = 1;
            break;
        case 'm':
            mapfilename = optarg;
            break;
        case 'n':
            npeople = atoi(optarg);
            break;
        case 'c':
            chrom = optarg;
            break;
        case 'o':
            outfilename = optarg;
            break;
        case 's':
            skipd = atoi(optarg);
            break;
        case 't':
            noutcomes = atoi(optarg);
            break;
        case 'g':
            ngpreds = atoi(optarg);
            break;
        case 'a':
            sep = optarg;
            break;
        case 'e':
            nohead = 1;
            break;
        case 'r':
            score = 1;
            break;
        case 'l':
            allcov = 1;
            break;
        case 'b':
            interaction = atoi(optarg);
            break;
        case 'k':
            interaction_excluded = atoi(optarg);
            break;
        case 'v':
            inverse_filename = optarg;
            break;
        case 'u':
            robust = 1;
            break;

        case '?':
            print_usage(program_name, 1);
        case -1:
            break;
        default:
            abort();
        } // end of switch
    } while (next_option != -1);
} // end of function

bool cmdvars::isIsInteractionExcluded() const
{
    return is_interaction_excluded;
}

void cmdvars::printinfo()
{
    fprintf(stdout,
            "%s v. %s (C) Yurii Aulchenko, Lennart C. Karssen, Maksim Struchalin, EMCR\n\n",
            PACKAGE, PACKAGE_VERSION);
#if EIGEN
    fprintf(stdout, "Using EIGEN for matrix operations\n");
#endif

    if (neco[0] != 1 || neco[1] != 1 || neco[2] != 1)
    {
        print_usage(program_name, 1);
    }

    if (score)
    {
        cout << "option --score suppressed from v 0.1-6\n";
        exit(1);
    }

    str_genfilename = genfilename;
    if (str_genfilename.find(FILEVECTOR_INDEX_FILE_SUFFIX) != string::npos
            || str_genfilename.find(FILEVECTOR_DATA_FILE_SUFFIX)
                    != string::npos)
        isFVF = 1;

    fprintf(stdout, "Options in effect:\n");
    fprintf(stdout, "\t --pheno   = %s\n", phefilename);
    fprintf(stdout, "\t --info    = %s\n", mlinfofilename);
    fprintf(stdout, "\t --dose    = %s\n", genfilename);
    if (isFVF)
        fprintf(stdout, "\t             (using FVF data)\n");
    fprintf(stdout, "\t --ntraits = %d\n", noutcomes);
    fprintf(stdout, "\t --ngpreds = %d\n", ngpreds);
    fprintf(stdout, "\t --interaction = %d\n", interaction);
    fprintf(stdout, "\t --interaction_only = %d\n", interaction_excluded);

    if (inverse_filename != NULL)
        fprintf(stdout, "\t --mmscore = %s\n", inverse_filename);
    else
        fprintf(stdout, "\t --mmscore     = not in output\n");
//  fprintf(stdout,"\t --mmscore = %s\n",inverse_filename);

    if (mapfilename != NULL)
        fprintf(stdout, "\t --map     = %s\n", mapfilename);
    else
        fprintf(stdout, "\t --map     = not in output\n");
    if (npeople > 0)
        fprintf(stdout, "\t --nids    = %d\n", npeople);
    else
        fprintf(stdout, "\t --nids    = estimated from data\n");
    if (chrom != "-1")
        cout << "\t --chrom   = " << chrom << "\n";
    else
        cout << "\t --chrom   = not in output\n";
    if (outfilename != NULL)
        fprintf(stdout, "\t --out     = %s\n", outfilename);
    else
        fprintf(stdout, "\t --out     = regression.out.txt\n");
    fprintf(stdout, "\t --skipd   = %d\n", skipd);
    cout << "\t --separat = \"" << sep << "\"\n";
    if (score)
        fprintf(stdout, "\t --score   = ON\n");
    else
        fprintf(stdout, "\t --score   = OFF\n");
    if (nohead)
        fprintf(stdout, "\t --nohead  = ON\n");
    else
        fprintf(stdout, "\t --nohead  = OFF\n");
    if (allcov)
        fprintf(stdout, "\t --allcov  = ON\n");
    else
        fprintf(stdout, "\t --allcov  = OFF\n");
    if (robust)
        fprintf(stdout, "\t --robust  = ON\n");
    else
        fprintf(stdout, "\t --robust  = OFF\n");

    if (ngpreds != 1 && ngpreds != 2)
    {
        fprintf(stderr,
                "\n\n--ngpreds should be 1 for MLDOSE or 2 for MLPROB\n");
        exit(1);
    }

    if (interaction_excluded != 0)
    {
        interaction = interaction_excluded; //ups
        is_interaction_excluded = true;
    }
    if (outfilename == NULL)
    {
        outfilename = (char *) string("regression").c_str();
    }
#if COXPH
    if (score)
    {
        fprintf(stderr, "\n\nOption --score is implemented for linear and logistic models only\n");
        exit(1);
    }

    if (inverse_filename != NULL)
    {
        std::cerr << "ERROR: mmscore is forbidden for cox regression\n";
        exit(1);
    }
    if (robust)
    {
        std::cerr << "ERROR: robust standard errors not implemented for Cox regression\n";
        exit(1);
    }
#endif
}
