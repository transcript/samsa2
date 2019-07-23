#!/usr/bin/env Python
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
# db_results_swapper.py
# Created 10/27/2017, this version edited 10/27/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program takes a list of IDs from one database, such as the RefSeq
# database, and gets all IDs for those entries retrieved against another 
# database, such as Subsystems.
#
# Usage: 
#
# -I		infile			specifies the infile (a DIAMOND results file 
#								in m8 format)
# -A		annotation		specifies annotated results to search against 
#
# -O		outfile			specifies a name for the outfile (otherwise defaults 
#								to $name_hierarchy.tsv)
# -P 		partial			partial outfile; a list of all reads with their
#								hierarchy results (OPTIONAL)
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

t0 = time.clock()

# loading starting file
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

infile = open (infile_name, "r")

# setting up storage variables
id_list = []
line_counter = 0

# reading through the infile
for line in infile:
	line_counter += 1
	splitline = line.split("\t")
	if line_counter % 1000000 == 0:
		t99 = time.clock()
		print str(line_counter)[:-6] + "MM lines processed so far in " + str(t99-t0) + " seconds."
	
	id_list.append(splitline[0])

t1 = time.clock()

# results reporting
print "\nAnalysis of " + infile_name + " complete."
print "Number of IDs found: " + str(len(id_list))
print "Time elapsed: " + str(t1-t0) + " seconds."

infile.close()
	
# time to search for these in the annotated file
if "-A" in sys.argv:
	af_name = string_find("-A")
else:
	sys.exit( "No annotation file indicated; skipping annotation filtering step.")

# IO
af = open (af_name, "r")
if "-O" in sys.argv:
	outfile_name = string_find("-O")
else:
	outfile_name = infile_name + ".subsys_results"

outfile = open (outfile_name, "w")

print "\nStarting annotation filtering now."

t2 = time.clock()

# data storage
af_line_counter = 0
af_hit_counter = 0
af_dic = {}

for line in af:
	af_line_counter += 1
	splitline = line.split("\t", 1)
	af_dic[splitline[0]] = splitline[1]
	if af_line_counter % 1000000 == 0:
		t3 = time.clock()
		print str(af_line_counter)[:-6] + "MM lines processed so far in " + str(t3-t2) + " seconds."

print "Annotation file read in."
entry_count = 0
for entry in id_list:
	entry_count += 1
	if entry_count % 1000000 == 0:
		print str(entry_count)[:-6] + "MM entries searched."
	try:
		af_entry = af_dic[entry]
		outfile.write(entry + "\t" + af_entry)
		af_hit_counter += 1
	except KeyError:
		continue

print "Comparison complete: " + str(af_hit_counter) + " hits found in " + str(af_line_counter) + " total lines."
print "Results saved as " + outfile_name + "."

af.close()
outfile.close()
