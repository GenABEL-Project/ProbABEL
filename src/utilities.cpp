/*
 * utilities.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: mkooyman
 *
 *
 * Copyright (C) 2009--2014 Various members of the GenABEL team. See
 * the SVN commit logs for more details.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301, USA.
 *
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
