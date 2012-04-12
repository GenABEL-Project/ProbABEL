/*
 * utilities.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: mkooyman
 */
#include <string>
void report_error(const char * format, ...)
{
    va_list args;
    char buffer[256];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    printf("ERROR: %s\n", buffer);
    exit(EXIT_FAILURE);
}

