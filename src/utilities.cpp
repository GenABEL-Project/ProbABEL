/*
 * utilities.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: mkooyman
 */
#include <string>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <cstdlib>

void report_error(const char * format, ...)
{
    va_list args;
    char buffer[256];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    std::cerr << "ERROR: " << buffer << std::endl;
    exit(EXIT_FAILURE);
}
