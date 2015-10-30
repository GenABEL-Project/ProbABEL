#include <map>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

#ifndef _NOT_R_FILEVECTOR
#include <R.h>
#endif
#include "frutil.h"
#include "CastUtils.h"

unsigned short int UNSIGNED_SHORT_INT_NAN;
short int SHORT_INT_NAN;
unsigned int UNSIGNED_INT_NAN;
int INT_NAN;
char CHAR_NAN;
unsigned char UNSIGNED_CHAR_NAN;
char const* parseFormats[9];


int initConsts(){
    int i;
    unsigned int ui;
    sscanf("32767", "%hi", &SHORT_INT_NAN);
    sscanf("65535", "%hu", &UNSIGNED_SHORT_INT_NAN);
    sscanf("2147483647", "%i", &INT_NAN);
    sscanf("4294967295", "%u", &UNSIGNED_INT_NAN);
    sscanf("127", "%i", &i); CHAR_NAN = i;
    sscanf("255", "%u", &ui); UNSIGNED_CHAR_NAN = ui;

    parseFormats[UNSIGNED_SHORT_INT] = "%hu";
    parseFormats[SHORT_INT] = "%hd";
    parseFormats[UNSIGNED_INT] = "%u";
    parseFormats[INT] = "%d";
    parseFormats[FLOAT] = "%f";
    parseFormats[DOUBLE] = "%lf";
    parseFormats[SIGNED_CHAR] = "%i";
    parseFormats[UNSIGNED_CHAR] = "%i";
    return 0;
}


int dummy = initConsts();


void parseStringToArbType(string s, int destType, void *destData,
                          string nanString) {

    char const *format = parseFormats[destType];

    int result;
    // no proper format specifier exists for char
    if (destType == SIGNED_CHAR || destType == UNSIGNED_CHAR) {
        int i;
        result = sscanf(s.c_str(), format, &i);
        if (nanString == s || result !=1){
            setNan(destData, destType);
            return;
        } else {
            if (destType == SIGNED_CHAR) *(char*) destData = i;
            if (destType == UNSIGNED_CHAR) *(unsigned char*) destData = i;
        }

    } else {
        result = sscanf(s.c_str(), format, destData);
        if (nanString == s || result !=1){
            setNan(destData, destType);
            return;
        }
    }
}


unsigned short int dataTypeFromString(string type){
    if (type == "UNSIGNED_SHORT_INT") return 1;
    if (type == "SHORT_INT") return 2;
    if (type == "UNSIGNED_INT") return 3;
    if (type == "INT") return 4;
    if (type == "FLOAT") return 5;
    if (type == "DOUBLE") return 6;
    if (type == "CHAR") return 7;
    if (type == "UNSIGNED_CHAR") return 8;
    return 0;
}


string dataTypeToString(int type){
    if (type == 1) return "UNSIGNED_SHORT_INT";
    if (type == 2) return "SHORT_INT";
    if (type == 3) return "UNSIGNED_INT";
    if (type == 4) return "INT";
    if (type == 5) return "FLOAT";
    if (type == 6) return "DOUBLE";
    if (type == 7) return "CHAR";
    if (type == 8) return "UNSIGNED_CHAR";
    return 0;
}


string bufToString(short int dataType, char *data, string nanString){
    char ret[500];
    switch(dataType){
    case UNSIGNED_SHORT_INT:
        sprintf(ret, "%hu", *(unsigned short int*)data);
        break;
    case SHORT_INT:
        sprintf(ret, "%hd", *(short int*)data);
        break;
    case UNSIGNED_INT:
        sprintf(ret, "%u", *(unsigned int*)data);
        break;
    case INT:
        sprintf(ret, "%d", *(int*)data);
        break;
    case FLOAT:
        sprintf(ret, "%f", *(float*)data);
        break;
    case DOUBLE: // changed to "%f" from %lf [not ISO C++]
        sprintf(ret, "%f", *(double*)data);
        break;
    case SIGNED_CHAR: // changed to "%f" from %lf [not ISO C++]
        sprintf(ret, "%d", (int)*(char*)data);
        break;
    case UNSIGNED_CHAR: // changed to "%f" from %lf [not ISO C++]
        sprintf(ret, "%d", (int)*(unsigned char*)data);
        break;
    }
    if (checkNan(data,dataType)) {
        return nanString;
    }

    return string(ret);
}

void setNan(unsigned short int &i){setNan(&i, UNSIGNED_SHORT_INT);}
void setNan(short int &i){setNan(&i, SHORT_INT);}
void setNan(unsigned int &i){setNan(&i, UNSIGNED_INT);}
void setNan(int &i){setNan(&i, INT);}
void setNan(float &i){setNan(&i, FLOAT);}
void setNan(double &i){setNan(&i, DOUBLE);}
void setNan(char &i){setNan(&i, SIGNED_CHAR);}
void setNan(unsigned char &i){setNan(&i, UNSIGNED_CHAR);}

bool checkNan(unsigned short int i){return checkNan(&i, UNSIGNED_SHORT_INT);}
bool checkNan(short int i){return checkNan(&i, SHORT_INT);}
bool checkNan(unsigned int i){return checkNan(&i, UNSIGNED_INT);}
bool checkNan(int i){return checkNan(&i, INT);}
bool checkNan(float i){return checkNan(&i, FLOAT);}
bool checkNan(double i){return checkNan(&i, DOUBLE);}
bool checkNan(char i){return checkNan(&i, SIGNED_CHAR);}
bool checkNan(unsigned char i){return checkNan(&i, UNSIGNED_CHAR);}

void setNan(void *data, int dataType){
    double dZero = 0.;
    float fZero = 0.;

    switch (dataType) {
    case UNSIGNED_SHORT_INT:
        (*(unsigned short int*) data) = UNSIGNED_SHORT_INT_NAN;
        break;
    case SHORT_INT:
        (*(short int*) data) = SHORT_INT_NAN;
        break;
    case UNSIGNED_INT:
        (*(unsigned int*) data) = UNSIGNED_INT_NAN;
        break;
    case INT:
        (*(int*) data) = INT_NAN;
        break;
    case FLOAT:
        (*(float*) data) = fZero/fZero;
        break;
    case DOUBLE:
        (*(double*) data) = dZero/dZero;
        break;
    case SIGNED_CHAR:
        (*(char*) data) = CHAR_NAN;
        break;
    case UNSIGNED_CHAR:
        (*(unsigned char*) data) = UNSIGNED_CHAR_NAN;
        break;
    default:
        errorLog << "file contains data of unknown type " << dataType
                 << endl << errorExit;
    }
}


bool checkNan(void *data, int dataType){
    switch (dataType) {
    case UNSIGNED_SHORT_INT:
        return (*(unsigned short int*) data) == UNSIGNED_SHORT_INT_NAN;
    case SHORT_INT:
        return (*(short int*) data) == SHORT_INT_NAN;
    case UNSIGNED_INT:
        return (*(unsigned int*) data) == UNSIGNED_INT_NAN;
    case INT:
        return (*(int*) data) == INT_NAN;
    case FLOAT:
#ifdef _NOT_R_FILEVECTOR
        return std::isnan(*(float*) data);
#else
        return ISNAN(*(float*) data);
#endif
    case DOUBLE:
#ifdef _NOT_R_FILEVECTOR
        return std::isnan(*(double*)data);
#else
        return ISNAN(*(double*)data);
#endif
    case UNSIGNED_CHAR:
        return (*(unsigned char*) data) == UNSIGNED_CHAR_NAN;
    case SIGNED_CHAR:
        return (*(char*) data) == CHAR_NAN;
    default:
        errorLog << "file contains data of unknown type " << dataType
                 << endl << errorExit;
        return false;
    }
}


int getDataType(unsigned short int){return UNSIGNED_SHORT_INT;}
int getDataType(short int){return SHORT_INT;}
int getDataType(unsigned int){return UNSIGNED_INT;}
int getDataType(int){return INT;}
int getDataType(float){return FLOAT;}
int getDataType(double){return DOUBLE;}
int getDataType(char){return SIGNED_CHAR;}
int getDataType(unsigned char){return UNSIGNED_CHAR;}
