# This AWK script verifies the beta and sebeta values for the flipmaf
# option of ProbABEL.
#
# Author: L.C. Karssen
# Copyright (c) 2016, PolyOmica, Groningen, The Netherlands,
# www.polyomica.com
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.
#


BEGIN {
    if (ARGC != 3) {
        print "Error: Please specify two files on the command line",
            "" > "/dev/stderr"
        exit 1
    }

    # Define offset of beta and se_beta columns. Calculate from last
    # column backwards.
    betacol_offset   = 2               # so that betacol =
                                       # $(NF-betacol_offset)
    sebetacol_offset = 1
}

# $4 is the Freq1 column denoting the frequency of the A1 allele.
NR != 1 && $4 > 0.5 {
    # $1 is the Rsq column; only check variants that have an
    # imputation Rsq > 1e-16, which is the threshold we use in
    # ProbABEL as well.
    if ($7 > 1e-16) {
        snpname = $1
        if ( !(snpname in beta) ) {
            # When reading the first file: the SNP entry shouldn't
            # exist so we add them to the array.

            # I'm using $(NF-x) here so that this check also works on
            # 2df files, which have more beta_SNP columns and only the
            # last ones are interesting.
            beta[snpname]   = $(NF-betacol_offset)
            sebeta[snpname] = $(NF-sebetacol_offset)
        } else {
            # Here we end up when reading the second file. This checks
            # if the current beta and sebeta are equal to the values
            # read from the first file.
            if (beta[snpname] == $(NF-betacol_offset) && beta[snpname] == "nan" &&
                sebeta[snpname] == $(NF-sebetacol_offset) && sebeta[snpname] == "nan")
            {
                # Beta and sebeta are nan in both files
                print snpname, "\tOK"
            }
            else if (beta[snpname] != -$(NF-betacol_offset) ||
                     sebeta[snpname != $(NF-sebetacol_offset)])
            {
                # In the 2df case of rare alleles, one beta,sebeta
                # combo is 0,0 and in that case we need to select the
                # other beta,sebeta columns
                if ($(NF-betacol_offset) == 0 && $(NF-sebetacol_offset) == 0)
                {
                    if (beta[snpname] == -$(NF-betacol_offset-2) &&
                        sebeta[snpname] == $(NF-sebetacol_offset-1))
                    {
                        print snpname, "\tOK"
                    }
                }
                else
                {
                print "Discordant results:", $0, "vs.",
                    beta[snpname], "+/-", sebeta[snpname] > "/dev/stderr"
                }
            }
            else
            {
                print snpname, "\tOK"
            }
        }
    }
}

END {
}
