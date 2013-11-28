/*
 * main_functions_dump.h
 *
 *  Created on: Nov 27, 2013
 *      Author: mkooyman
 */

#ifndef MAIN_FUNCTIONS_DUMP_H_
#define MAIN_FUNCTIONS_DUMP_H_

#include "phedata.h"
#include "command_line_settings.h"
#include "maskedmatrix.h"

void update_progress_to_cmd_line(int, int);
void loadInvSigma(cmdvars& input_var, phedata& phd,
		masked_matrix& invvarmatrix);
int create_phenotype(phedata& phd, cmdvars& input_var);
void create_start_of_header(std::vector<std::ofstream*>& outfile,
		cmdvars& input_var, phedata& phd);
void write_mlinfo(const std::vector<std::ofstream*>& outfile, unsigned int file,
		const mlinfo& mli, int csnp, const cmdvars& input_var, int gcount,
		double freq);
void open_files_for_output(std::vector<std::ofstream*>& outfile,
		std::string& outfilename_str);
void create_header(std::vector<std::ofstream*>& outfile, cmdvars& input_var,
		phedata& phd, int& interaction_cox);
int get_start_position(const cmdvars& input_var, int model,
        int number_of_rows_or_columns);
#endif /* MAIN_FUNCTIONS_DUMP_H_ */
