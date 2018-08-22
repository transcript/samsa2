#!/usr/lib/python2.7

##########################################################################
#
# Copyright (C) 2015-2016 Sam Westreich
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation;
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
##########################################################################
#
# subsys_reducer.py
# Created 3/01/2017, this version created 5/22/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program parses through the results file from a DIAMOND annotation run
# against the Subsystems database, in BLAST m8 format, to get the results
# into a more compressed format for statistical analysis.
#
# Usage:
#
# -I		infile			specifies the infile (a DIAMOND results file
#								in m8 format)
# -O		outfile			optional; default outfile name is infile.reduced
#
##########################################################################

# imports
import sys

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem

# infile
try:
	infile = open(string_find("-I"), "r")
except IndexError:
	sys.exit("Unable to identify infile.  Is it specified with '-I' flag?")

# specifying variables
percent_dic = {}
count_dic = {}
line_counter = 0
unique_l4 = {}
upper_lvls = {}

# parsing
for line in infile:
	line_counter += 1
	splitline = line.split("\t", 2)
	split_levels = splitline[2].split("\t")

	if split_levels[0] not in unique_l4.keys():
		unique_l4[split_levels[0]] = float(splitline[1])
		upper_lvls[split_levels[0]] = splitline[2]
	else:
		if float(splitline[1]) >= unique_l4[split_levels[0]]:
			unique_l4[split_levels[0]] = float(splitline[1])
			upper_lvls[split_levels[0]] = splitline[2]
		else:
			continue

infile.close()

# at this point, the infile is reopened and parsed a second time, this time
# against the previously assembled parsed version.

infile = open(string_find("-I"), "r")

# parsing, round two
for line in infile:
	splitline = line.split("\t")

	try:
		percent_dic[splitline[2]] += float(splitline[0])
	except KeyError:
		percent_dic[splitline[2]] = float(splitline[0])

	try:
		count_dic[splitline[2]] += float(splitline[1])
	except KeyError:
		count_dic[splitline[2]] = float(splitline[1])

# creating the outfile
if "-O" in sys.argv:
	outfile_name = string_find("-O")
else:
	outfile_name = string_find("-I") + ".reduced"
outfile = open(outfile_name, "w")

# writing out to outfile
for k, v in sorted(count_dic.items(), key=lambda kv: -kv[1]):
	outfile.write(str(percent_dic[k]) + "\t" + str(count_dic[k]) + "\t" + upper_lvls[k])

outfile.close()

# reporting final statistics
print ("Total number of lines analyzed:\t" + str(line_counter))
print ("Total number of final entries:\t" + str(len(percent_dic.keys())))
print ("Results saved in file:\t\t" + outfile_name)
