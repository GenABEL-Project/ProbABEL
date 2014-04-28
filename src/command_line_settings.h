/*
 * command_line_settings.h
 *
 *  Created on: Apr 2; int 2012
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


#ifndef COMMAND_LINE_SETTINGS_H_
#define COMMAND_LINE_SETTINGS_H_
#include <string>

using std::string;

class cmdvars
{
 private:
    char * program_name;

    char *phefilename;
    char *mlinfofilename;
    char *genfilename;
    char *mapfilename;
    string outfilename;
    char *inverse_filename;

    string str_genfilename;

    int nohead;
    int score;
    int npeople;
    int ngpreds;
    int interaction;
    int interaction_excluded;
    bool is_interaction_excluded;
    int robust;
    string chrom;
    string sep;
    int neco[3];                /* Necessary command line options */
    bool iscox;
    int isFVF;
    int noutcomes;
    int skipd;
    int allcov;

 public:
    cmdvars()
    {
        program_name = NULL;

        std::fill_n(neco, 3, 0);
        phefilename      = NULL;
        mlinfofilename   = NULL;
        genfilename      = NULL;
        mapfilename      = NULL;
        outfilename      = string("");
        inverse_filename = NULL;

        sep     = " ";
        nohead  = 0;
        score   = 0;
        npeople = -1;
        ngpreds = 1;
        interaction = 0;
        interaction_excluded = 0;
        is_interaction_excluded = false; //Oh Holy Matrix, forgive me for this!
        robust = 0;
        chrom = "-1";
        str_genfilename = "";
        isFVF  = 0;
        skipd  = 2;
        allcov = 0;
#if COXPH
        noutcomes = 2;
        iscox     = true;
#else
        noutcomes = 1;
        iscox     = false;
#endif
    }

    void set_variables(int, char *[]);
    char* getPhefilename() const;
    int getAllcov() const;
    string getChrom() const;
    char* getGenfilename() const;
    int getInteraction() const;
    char* getInverseFilename() const;
    bool isIscox() const;
    int getIsFvf() const;
    char* getMapfilename() const;
    char* getMlinfofilename() const;

    int getNgpreds() const;
    int getNohead() const;
    int getNoutcomes() const;
    int getNpeople() const;
    string getOutfilename() const;
//TODO(unknown) This function is not used. Remove in near future
//    char* getProgramName() const;
    int getRobust() const;
    int getScore() const;
    string getSep() const;
    int getSkipd() const;
    string getStrGenfilename() const;

    void printinfo();
    bool isIsInteractionExcluded() const;
};

#endif /* COMMAND_LINE_SETTINGS_H_ */
