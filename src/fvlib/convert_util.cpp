#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace std;

#include "FileVector.h"
#include "Transposer.h"

#include "convert_util.h"

#define REPORT_EVERY 10000
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

unsigned long calcNumLines(string fileName){
    ifstream file(fileName.c_str());
    string line;
    unsigned long numlines = 0;
    while(getline(file, line)){
        numlines++;
    }
    return numlines;
}


void messageOnOff(int bit)
{
    msg << (bit?"ON":"OFF");
}


unsigned long calcNumWordsInFirstLine(string fileName){
    ifstream file(fileName.c_str());
    vector<string> words;
    string line;
    getline(file, line);
    tokenize(line, words);
    return words.size();
}


void text2fvf(string program_name, string infilename, string outfilename,
              string rownamesfilename, string colnamesfilename,
              int rncol, int cnrow,
              unsigned long skiprows, unsigned long skipcols,
              int bTranspose, int Rmatrix, unsigned short type,
              bool quiet, string nanString) {

    if (Rmatrix) {
        skipcols = skiprows = 1;
        rncol = cnrow = 1;
    }

    if (!quiet) {
        msg << "Options in effect: \n";

        msg << "\t --infile    = " << infilename.c_str() << endl;
        msg << "\t --outfile   = " << outfilename.c_str()<< endl;

        msg << "\t --skiprows  = ";
        if (skiprows) msg << skiprows << endl; else msg << "OFF\n";

        msg << "\t --skipcols  = ";
        if (skipcols) msg << skipcols << endl; else msg << "OFF\n";

        msg << "\t --cnrow     = ";
        if (cnrow || colnamesfilename != "") {
            if (colnamesfilename == "")
                msg << "ON, using line " << cnrow << " of '"
                    << infilename << "'\n";
            else
                msg << "ON, using data from file '"
                    << colnamesfilename<< "'\n";
        } else {
            msg << "OFF\n";
        }

        msg << "\t --rncol     = ";

        if (rncol || rownamesfilename != "") {
            if (rownamesfilename == "")
                msg << "ON, using column " << rncol << " of '"
                    << infilename << "'\n";
            else
                msg << "ON, using data from file '"
                    << rownamesfilename << "'\n";
        } else {
            msg << "OFF\n";
        }

        msg << "\t --transpose = ";
        messageOnOff(bTranspose);
        msg << "\n";
        msg << "\t --Rmatrix   = ";
        messageOnOff(Rmatrix);
        msg << endl;
        msg << "\t --nanString = " << nanString;
        msg << endl;
    }

    const string TMP_SUFFIX = "_fvtmp";
    string realOutFilename = outfilename;

    // if transpose if OFF, do transpose (sic) and write to temp file
    if (!bTranspose){
        outfilename = outfilename + TMP_SUFFIX;
    }

    // check that it is mentioned how many columns to skip when reading row names
    if (rownamesfilename=="" && (rncol > 0) && !skipcols) {
        errorLog << "\n\nPlease tell how many columns to skip when you "
                 << "provide the column containing row names!"
                 << endl << endl << errorExit;
    }
    if (skipcols && (rncol > ((int) skipcols) ) ) {
        errorLog << "rncol > skipcols" << errorExit;
    }
    if (colnamesfilename=="" && (cnrow>0) && !skiprows) {
        errorLog << "\n\nPlease tell how many rows to skip when you "
                 << "provide the row containing column names!"
                 << endl << endl << errorExit;
    }
    if (skiprows && (cnrow > ((int) skiprows) ) ) {
        errorLog << "cnrow > skiprows" << errorExit;
    }

    ifstream infile(infilename.c_str());
    if (!infile) {
        errorLog << "Can not open file '" << infilename
                 << "' for reading\n\n" << errorExit;
    }

    vector<string> extColNames;

    // column names specified in a separate file; check the number of columns
    unsigned long words_in_colnamesfile = 0;

    // reading column names from file
    if (colnamesfilename != "")
    {
        msg << "Reading columns from "  << colnamesfilename << ": ";
        ifstream colnamesfile(colnamesfilename.c_str());
        string tmpstr;
        if (!colnamesfile) {
            errorLog << "Can not open column names file '"
                     << colnamesfilename << "'\n\n" << errorExit;
        }
        while (colnamesfile >> tmpstr) {
            words_in_colnamesfile++;
            extColNames.push_back(tmpstr);
        }
        colnamesfile.close();
        if (!quiet) {
            msg << "Number of names in column-names file '"
                << colnamesfilename << "' is ";
            msg << words_in_colnamesfile<< "\n\n";
        }
    }

    vector<string> extRowNames;

    // row names specified in a separate file; check the number of columns
    unsigned long words_in_rownamesfile = 0;

    // reading row names from file
    if (rownamesfilename != "") {
        ifstream rownamesfile(rownamesfilename.c_str());
        string tmpstr;
        if (!rownamesfile) {
            errorLog << "Can not open row names file '"
                     << rownamesfilename << "'\n\n";
        }
        while (rownamesfile >> tmpstr) {
            words_in_rownamesfile++;
            extRowNames.push_back(tmpstr);
        }

        rownamesfile.close();

        if (!quiet){
            msg << "number of names in row-names file '" << rownamesfilename;
            msg << "' is " << words_in_rownamesfile << "\n\n" << endl;
        }

    }

    ifstream srcFile(infilename.c_str());

    string firstLine;

    vector<string> firstLineWords;
    tokenize(firstLine, firstLineWords);

    unsigned long numLines = calcNumLines(infilename);
    unsigned long numWords = calcNumWordsInFirstLine(infilename);

    msg << "Number of lines in source file is " << numLines << endl;
    msg << "Number of words in source file is " << numWords << endl;

    msg << "skiprows = " << skiprows << endl;
    msg << "cnrow = "    << cnrow    << endl;
    msg << "skipcols = " << skipcols << endl;
    msg << "rncol = "    << rncol    << endl;
    msg << "Rmatrix = "  << Rmatrix  << endl;
    msg << "numWords = " << numWords << endl;

    long int numRows = numLines - skiprows;
    long int numColumns=123;

    bool colNamesFilePresents = (colnamesfilename!="");
    bool rowNamesFilePresents = (rownamesfilename!="");

    numColumns = numWords - skipcols;
    if (Rmatrix && !colNamesFilePresents)
        numColumns = numWords - skipcols + 1;

    msg << "Creating file with numRows = " << numRows << endl;
    msg << "Creating file with numColumns = " << numColumns << endl;

    initializeEmptyFile(outfilename, numRows, numColumns, type, true );
    FileVector *out = new FileVector(outfilename, 1);

    string line;
    unsigned long rowCnt = 1;
    unsigned long lineCnt = 1;

    short element_size = out->getElementSize();
    char* ArbTypeData =
        new (nothrow) char [(out->getNumObservations()) * element_size];

    while(getline(srcFile, line)) {

        // is this a column name line?
        if (lineCnt == ((unsigned long) cnrow) && !colNamesFilePresents)
        {
            unsigned long i;
            vector<string> lineWords;
            tokenize(line, lineWords);
            // ignoring R-matrix flag for some reason
            for (i = skipcols-Rmatrix; i < lineWords.size(); i++) {
                extColNames.push_back(lineWords[i]);
            }
            lineCnt++;
            continue;
        }

        if (lineCnt <= skiprows)
        {
            lineCnt++;
            continue;
        }

        unsigned long colCnt = 0;
        unsigned long i = 0;
        char const *SEPARATORS = " \t";
        char *lineBuf = strdup(line.c_str());
        char *word = strtok(lineBuf, SEPARATORS);
        while(word!=0) {

            if (i == ( (unsigned long) (rncol - 1) ) && !rowNamesFilePresents) {
                extRowNames.push_back(word);
                word = strtok(0, SEPARATORS);
                i++;
                continue;
            }

            if (i < skipcols) {
                word = strtok(0, SEPARATORS);
                i++;
                continue;
            }

            parseStringToArbType(word, type,
                                 &ArbTypeData[colCnt*element_size], nanString);
            colCnt++;
            word = strtok(0, SEPARATORS);
            i++;
        }
        free(lineBuf);

        out->writeVariable(rowCnt-1, (char*)ArbTypeData);
        rowCnt++;
        lineCnt++;
    }

    delete [] ArbTypeData;

    if (!colNamesFilePresents && cnrow < 0)
    {
        for (long int i = 1; i <= numColumns; i++) {
            char * tmpstr = NULL;
            sprintf(tmpstr,"%lu",i);
            extColNames.push_back(tmpstr);
        }
    }
    if (!rowNamesFilePresents && cnrow < 0)
    {
        for (long int i = 1; i <= numRows; i++) {
            char * tmpstr = NULL;
            sprintf(tmpstr,"%lu",i);
            extRowNames.push_back(tmpstr);
        }
    }
    unsigned long i;
    for (i = 0; i < extColNames.size(); i++){
        out->writeObservationName(i, extColNames[i]);
    }

    for (i = 0; i < extRowNames.size(); i++){
        out->writeVariableName(i, extRowNames[i]);
    }

    delete out;

    if(!bTranspose)
    {
        msg << "Transposing " << outfilename << " => "
            << realOutFilename << "." << endl;
        Transposer tr;

        tr.process(outfilename, realOutFilename, true);
    }
    msg << "text2fvf finished." << endl;
}
