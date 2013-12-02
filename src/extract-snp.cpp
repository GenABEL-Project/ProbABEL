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
         << "file in filevector format (.fvi/.fvd)." << endl;
    cout << endl;
    cout << "Usage: " << program_name << " --file <fv file> --snp <snpname>"
         << endl;
    cout << "   or: " << program_name << " -f <fv file> -s <snpname>" << endl;
    cout << endl;
    cout << "Additional options:" << endl;
    cout << "\t--id (-i) <idname>: "
         << "only print the dosage for the given SNP and ID" << endl;
    cout << "\t--debug (-d): show debugging output" << endl;
    cout << "\t--help (-h): show this information" << endl;
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        info(argv[0]);
        exit(3);
    }

    int next_option;
    const char * const short_options = "df:hi:s:";
    const struct option long_options [] =
        {
            {"debug",  0, NULL, 'd'},
            {"file",  1, NULL, 'f'},
            {"help",  0, NULL, 'h'},
            {"id",   1, NULL, 'i'},
            {"snp",  1, NULL, 's'},
            {NULL  ,   0, NULL, 0  }
        };
    char *program_name = argv[0];
    bool debug = false;
    string inputFileName = "";
    string idname = "";
    string snpname = "";
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
        case 's':
            snpname = string(optarg);
            break;
        case '?': break;
        case -1 : break;
        }
    } while (next_option != -1);

    if (snpname.compare("") == 0)
    {
        cerr << "Error: No SNP name given" << endl << endl;
        info(program_name);
        exit(2);
    }

    unsigned long int snprow = 0;
    bool snpfound = false;
    bool idfound  = false;

    if (debug) cout << "Input file is '" << inputFileName << "'." << endl;

    FileVector fv(inputFileName, 64, true);

    if (debug)
    {
        for (unsigned long int col = 0; col < fv.getNumObservations(); col++)
        {
            cout << fv.readObservationName(col).name << " ";
        }

        cout << endl;
        cout << "----------------------" << endl;
    }

    // Look at the SNPs (rows) first
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
        cout << "----------------------" << endl;

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
