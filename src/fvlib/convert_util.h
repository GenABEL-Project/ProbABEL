#ifndef __CONVERT_UTIL__
#define __CONVERT_UTIL__

using namespace std;

unsigned long calcNumLines(string fileName);
void parseStringToArbType(string s, int destType, void *destData);
void text2fvf(string program_name, string infilename, string outfilename,
              string rownamesfilename, string colnamesfilename,
              int rncol, int cnrow,
              unsigned long skiprows, unsigned long skipcols,
              int bTranspose, int Rmatrix, unsigned short type,
              bool quiet, string nanString);

#endif
