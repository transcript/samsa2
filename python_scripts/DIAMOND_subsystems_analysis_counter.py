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
# DIAMOND_subsystems_analysis_counter.py
# Created 2/01/2017, this version edited 3/20/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program parses through the results file from a DIAMOND annotation run
# (in BLAST m8 format) to get the results into something more compressed
# and readable, against the SUBSYSTEMS database.
#
# Usage:
#
# -I		infile			specifies the infile (a DIAMOND results file
#								in m8 format)
# -O		outfile			specifies a name for the outfile (otherwise defaults
#								to $name_hierarchy.tsv)
# -P 		partial			partial outfile; a list of all reads with their
#								hierarchy results (OPTIONAL)
#
##########################################################################

# imports
import operator, sys, gzip, re

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			return next_elem

# loading starting file
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

infile = open (infile_name, "r")

# setting up databases
hit_count_db = {}
unique_seq_db = {}
read_id_db = {}
line_counter = 0

print ("\nNow reading through the SEED m8 results infile.")

# reading through the infile
for line in infile:
	line_counter += 1
	splitline = line.split("\t")

	unique_seq_db[splitline[0]] = 1

	if "-P" in sys.argv:
		read_id_db[splitline[0]] = splitline[1]

	try:
		hit_count_db[splitline[1]] += 1
	except KeyError:
		hit_count_db[splitline[1]] = 1
		continue



# results reporting
print ("\nAnalysis of " + infile_name + " complete.")
print ("Number of total lines: " + str(line_counter))
print ("Number of unique sequences: " + str(len(unique_seq_db)))

infile.close()

# IO
# try:
db = open (infile_name, "r", encoding='utf-8', errors='ignore')
# except TypeError:
# 	# error catching for Python 2
# 	db = open (db_name, "r")

if "-P" in sys.argv:
	partial_outfile_name = string_find("-P")
	partial_outfile = open(partial_outfile_name, "w")

print ("\nStarting database analysis now.")

# building a dictionary of the reference database
db_hier_dictionary = {}
db_line_counter = 0
db_error_counter = 0

for line in db:
	line = '\t'.join(line.split("\t")[12:])

	db_line_counter += 1
	splitline = line.split("\t", 1)

	# ID, the hit returned in DIAMOND results
	db_id = str(splitline[0]) # [1:]

	# name and functional description
	if "NO HIERARCHY" in splitline[1]:
		db_hier = "NO HIERARCHY"
	else:
		hier_split = splitline[1].split("\t")
		if hier_split[3].strip() != "":
			db_hier = hier_split[0] + "\t" + hier_split[1] + "\t" + hier_split[2] + "\t" + hier_split[3]
		else:
			db_hier = hier_split[0] + "\t" + hier_split[1] + "\t\t" + hier_split[2] + "\t" + hier_split[3]

	# add to dictionaries
	db_hier_dictionary[db_id] = db_hier

print ("\nSuccess!")
print ("Number of lines: " + str(db_line_counter))
print ("Number of errors: " + str(db_error_counter))

# printing out the partial outfile
if "-P" in sys.argv:
	for entry in read_id_db.keys():
		partial_outfile.write(entry + "\t" + read_id_db[entry] + "\t" + db_hier_dictionary[read_id_db[entry]] + "\n")

# condensing down the identical matches
condensed_hit_db = {}

for entry in hit_count_db.keys():
	org = db_hier_dictionary[entry]
	if org in condensed_hit_db.keys():
		condensed_hit_db[org] += hit_count_db[entry]
	else:
		condensed_hit_db[org] = hit_count_db[entry]

# dictionary output and summary
print ("\nDictionary database assembled.")
print ("Number of errors: " + str(db_error_counter))

print ("\nTop ten hierarchy matches:")
for k, v in sorted(condensed_hit_db.items(), key=lambda kv: -kv[1])[:10]:
	try:
		print (str(v) + "\t" + k )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue

# creating the outfiles
if "-O" in sys.argv:
	outfile_name = string_find("-O")
else:
	outfile_name = infile_name[:-4] + ".hierarchy"

outfile = open (outfile_name, "w")

# writing the output
error_counter = 0
for k, v in sorted(condensed_hit_db.items(), key=lambda kv: -kv[1]):
	try:
		q = v * 100 / float(line_counter)
		outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
	except KeyError:
		outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue

print ("\nAnnotations saved to file: '" + outfile_name + "'.")
print ("Number of errors: " + str(error_counter))

db.close()
outfile.close()
