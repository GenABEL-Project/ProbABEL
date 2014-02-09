/*
 * main_functions_dump.h
 *
 *  Created on: Nov 27, 2013
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


#ifndef MAIN_FUNCTIONS_DUMP_H_
#define MAIN_FUNCTIONS_DUMP_H_

#include "phedata.h"
#include "command_line_settings.h"
#include "maskedmatrix.h"

void update_progress_to_cmd_line(const int csnp, const int nsnps);

void loadInvSigma(const cmdvars& input_var, phedata& phd,
                  masked_matrix& invvarmatrix);

int create_phenotype(phedata& phd, const cmdvars& input_var);

void create_start_of_header(std::vector<std::ofstream*>& outfile,
                            cmdvars& input_var, phedata& phd);

void write_mlinfo(const std::vector<std::ofstream*>& outfile,
                  const unsigned int file, const mlinfo& mli,
                  const int csnp, const cmdvars& input_var,
                  const int gcount, const double freq);

void open_files_for_output(std::vector<std::ofstream*>& outfile,
                           const std::string& outfilename_str);

void create_header(std::vector<std::ofstream*>& outfile,
                   cmdvars& input_var, phedata& phd, int& interaction_cox);

int get_start_position(const cmdvars& input_var, const int model,
        const int number_of_rows_or_columns);
#endif /* MAIN_FUNCTIONS_DUMP_H_ */
