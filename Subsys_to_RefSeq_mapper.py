#!/usr/lib/python2.7
##########################################################################
#
# Copyright (C) 2015-2017 Sam Westreich
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
# Subsystems_specific_read_retriever.py
# Created 6/03/2017, this version updated 6/08/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# Purpose: Because all reads processed in the SAMSA database maintain their
# original read ID, it is possible to get results across databases!  For
# example, we can use this program to obtain all reads matching "fucosidase"
# from the Subsystems results and determine their organism of origin in the
# RefSeq results.
#
# This program takes a search term, the Subsystems ".receipt" file, and the 
# RefSeq results file as its inputs.  It checks the Subsystems ".receipt"
# file for the search term, finds the corresponding annotations in the RefSeq
# results, and returns those sequences.

# USAGE:
# -S		SEED results			The SEED Subsystems annotation for that sample.
# -I		Original input			The RefSeq results for that sample.
# -O		Output name				Name of saved retrieved sequences file
# -T		search Term				Term to be searched for in results for filtering
# 
##########################################################################

# imports
import operator, sys, time, gzip, re

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem

if "-I" not in sys.argv:
	sys.exit("Specify the file of reads to be retrieved with the -I flag.")
if "-S" not in sys.argv:
	sys.exit("Specify the SEED Subsystems annotation results with the -S flag.")
if "-O" not in sys.argv:
	sys.exit("Specify an output name with the -O flag.")
if "-T" in sys.argv:
	search_flag = True
	search_term = string_find("-T")
else:
	search_flag = False

# variables
read_id_list = []
hit_flag = False
i = 0
noMatch = 0
original_db = {}
line_counter = 0

# getting the read IDs from the original reads
sr_file = open(string_find("-I"), "r")
for line in sr_file:
	splitline = line.split("\t")
	if search_flag == True:
		if search_term in line:
			read_id_list.append(splitline[0])
	else:
		read_id_list.append(splitline[0])
	if i < 10:
		i += 1
#		print splitline[0]
	
	line_counter += 1
	if line_counter % 100000 == 0:
		print str(line_counter)[:-3] + "k lines analyzed so far."

sr_file.close()

print "\nRead IDs scanned: " + str(len(read_id_list)) + " IDs found.\n"

# reading in the original file
outfile = open(string_find("-O"), "w")
infile = open(string_find("-S"), "r")
i = 0
line_counter = 0

for line in infile:
	splitline = line.split("\t")
	if i < 10:
		i += 1
#		print splitline[0]
	original_db[splitline[0]] = line

	line_counter += 1
	if line_counter % 1000000 == 0:
		print str(line_counter)[:-6] + "M lines analyzed so far."

# matching and printing to output
in_common = set(original_db.keys()).intersection(read_id_list)
print "\n" + str(len(in_common)) + " reads in common.\nSaving..."

for entry in in_common:
	outfile.write(original_db[entry])

infile.close()
outfile.close()

print "Done scanning the input file; results are saved as " + string_find("-O")
print "Number of read IDs not found in original file: " + str(len(read_id_list)-len(in_common))