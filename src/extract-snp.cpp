/*******************************************************************************
 *
 * This program extracts a given row (dosages for all individuals for a given
 * SNP) from filevector files.
 *
 * (C) 2012 L.C. Karssen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <stdio.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#if HAVE_CONFIG_H
#include "config.h"
#endif

using std::cout;
using std::cerr;
using std::endl;

#include "fvlib/Logger.h"
#include "fvlib/FileVector.h"
#include "fvlib/CastUtils.h"


void info(char *program_name)
{
    cout << program_name
         << " extracts the dosage for a SNP for one or more individuals from a"
         << " file in filevector format (.fvi/.fvd)." << endl;
    cout << endl;
    cout << "Usage: " << program_name << " --file <fv file> --snp <snpname>"
         << endl;
    cout << "   or: " << program_name << " -f <fv file> -s <snpname>" << endl;
    cout << endl;
    cout << "Additional options:" << endl;
    cout << "\t--listsnps (-L); list all SNP names from the fv file"
         << " (the program exist after this)"
         << endl;
    cout << "\t--listids (-I); list all ID names from the fv file"
         << " (the program exist after this)"
         << endl;
    cout << "\t--id (-i) <idname>: "
         << "only print the dosage for the given SNP and ID" << endl;
    cout << "\t--debug (-d): show debugging output" << endl;
    cout << "\t--help (-h): show this information" << endl;
    cout << endl;
    cout << program_name << " is part of " << PACKAGE
         << " v" << PACKAGE_VERSION << endl;
    cout << "(C) 2015 Lennart C. Karssen, PolyΩmica, NL" << endl;
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        info(argv[0]);
        exit(3);
    }

    int next_option;
    const char * const short_options = "df:hi:ILs:";
    const struct option long_options [] =
        {
            {"debug",    0, NULL, 'd'},
            {"file",     1, NULL, 'f'},
            {"help",     0, NULL, 'h'},
            {"id",       1, NULL, 'i'},
            {"listids",  0, NULL, 'I'},
            {"listsnps", 0, NULL, 'L'},
            {"snp",      1, NULL, 's'},
            {NULL,       0, NULL, 0  }
        };

    char *program_name = argv[0];
    bool debug    = false;
    bool listSNPs = false;
    bool listIDs  = false;
    string inputFileName = "";
    string idname        = "";
    string snpname       = "";
    do
    {
        next_option = getopt_long(argc, argv,
                                  short_options, long_options, NULL);
        switch (next_option)
        {
        case 'h':
            info(program_name);
            exit(0);
        case 'd':
            debug = true;
            break;
        case 'f':
            inputFileName = string(optarg);
            break;
        case 'i':
            idname = string(optarg);
            break;
        case 'I':
            listIDs = true;
            break;
        case 'L':
            listSNPs = true;
            break;
        case 's':
            snpname = string(optarg);
            break;
        case '?': break;
        case -1 : break;
        }
    } while (next_option != -1);

    if (inputFileName.compare("") == 0)
    {
        cerr << "Error: No filevector file specified.\n";
        info(program_name);
        exit(4);
    }

    if (debug) cout << "Input file is '" << inputFileName << "'." << endl;


    FileVector fv(inputFileName, 64, true);

    if (listSNPs)
    {
        // Only print SNPs to stdout.
        for (unsigned long int row = 0; row < fv.getNumVariables(); row++)
        {
            string current_snpname = fv.readVariableName(row).name;
            cout << current_snpname << " ";
        }
        cout << endl;
        exit(0);
    }

    if (listIDs || debug)
    {
        if(debug)
        {
            cout << "---------- List of IDs ----------" << endl;
        }
        // Only print IDs to stdout.
        for (unsigned long int col = 0; col < fv.getNumObservations(); col++)
        {
            cout << fv.readObservationName(col).name << " ";
        }

        cout << endl;

        if(debug)
        {
            cout << "---------- End of IDs ----------" << endl;
        }

        if(!debug) exit(0);
    }


    if (snpname.compare("") == 0)
    {
        cerr << "Error: No SNP name given" << endl << endl;
        info(program_name);
        exit(2);
    }

    unsigned long int snprow = 0;
    bool snpfound = false;
    bool idfound  = false;


    // Look at the SNPs (rows) first
    if(debug)
    {
        cout << "---------- List of SNP names ----------" << endl;
    }
    for (unsigned long int row = 0; row < fv.getNumVariables(); row++)
    {
        string current_snpname = fv.readVariableName(row).name;
        if (current_snpname.compare(snpname) == 0)
        {
            snpfound = true;
            snprow = row;
            if (debug)
            {
                cout << "*" << current_snpname << "* ";
            }
        } else {
            if (debug)
            {
                cout << current_snpname << " ";
            }
        }
    }
    if (debug)
    {
        cout << endl;
        cout << "---------- End of SNP names ----------" << endl;

        cout << "N_obs = " << fv.getNumObservations() << "\telement size: "
             <<  fv.getElementSize() << "\tDataType: "
             << dataTypeToString(fv.getElementType()) << endl;
    }

    if (!snpfound)
    {
        cerr << "SNP name '" << snpname << "' not found in data file "
             << inputFileName << endl;
        exit(1);
    }

    char * data = new (nothrow) char[fv.getNumObservations() *
                                         fv.getElementSize()];
    if (!data)
    {
        cerr << "Cannot allocate memory for data vector. Exiting..." << endl;
        exit(2);
    }

    fv.readVariable(snprow, data);
    for (unsigned long int col = 0; col < fv.getNumObservations(); col++)
    {
        if (idname.compare("") != 0)
        {
            // An ID name has been set, only print the dosage for that ID
            // and SNP combination
            if (idname.compare(fv.readObservationName(col).name) == 0)
            {
                idfound = true;
                cout << fv.readObservationName(col).name << "\t"
                     << bufToString(fv.getElementType(),
                                    &data[col*fv.getElementSize()],
                                    string("NA"))
                     << endl;
            }
        } else {
            // Print the dosages for all IDs
            cout << fv.readObservationName(col).name << "\t"
                 << bufToString(fv.getElementType(),
                                &data[col*fv.getElementSize()],
                                string("NA"))
                 << endl;
        }
    }

    if (idfound == false && idname.compare("") != 0)
    {
        cerr << "Id '" << idname << "' not found in data file "
             << inputFileName << endl;
        exit(1);
    }

    delete [] data;
}
