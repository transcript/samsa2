#!/usr/bin/env Python
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
# DIAMOND_analysis_counter.py
# Created 8/16/2016, this version created 1/10/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program parses through the results file from a DIAMOND annotation run
# (in BLAST m8 format) to get the results into something more compressed
# and readable.
#
# Usage:
#
# -I		infile			specifies the infile (a DIAMOND results file
#								in m8 format)
# -O		organism		returns organism results
# -F		function		returns functional results
# -SO		specific org	creates a separate outfile for results that hit
#							a specific organism
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

# checking for an option (organism or function) to be specified
if "-O" not in sys.argv:
	if "-F" not in sys.argv:
		sys.exit("WARNING: need to specify either organism results (with -O flag in command) or functional results (with -F flag in command).")

# loading starting file
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

infile = open (infile_name, "r")

# setting up databases
RefSeq_hit_count_db = {}
unique_seq_db = {}
line_counter = 0

# reading through the infile - the DIAMOND results m8 format
print ("\nNow reading through the m8 results infile.")

for line in infile:
	line_counter += 1
	splitline = line.split("\t")
	# if line_counter % 1000000 == 0:
	# 	print (str(line_counter)[:-6] + "M lines processed so far in " + str(t99-t0) + " seconds.")

	unique_seq_db[splitline[0]] = 1

	try:
		RefSeq_hit_count_db[splitline[1]] += 1
	except KeyError:
		RefSeq_hit_count_db[splitline[1]] = 1
		continue


print ("\nAnalysis of " + infile_name + " complete.")
print ("Number of total lines: " + str(line_counter))
print ("Number of unique sequences: " + str(len(unique_seq_db)))

infile.close()

# time to search for these in the reference database
db = open (infile_name, "r")

print ("\nStarting database analysis now.")

# optional outfile of specific organism results
# if "-SO" in sys.argv:
# 	target_org = string_find("-SO")
# 	db_SO_dictionary = {}

# building a dictionary of the reference database
if "-F" in sys.argv:
	db_func_dictionary = {}
if "-O" in sys.argv:
	db_org_dictionary = {}
db_line_counter = 0
db_error_counter = 0

# Reusing db code for parsing diamond file with additional salltitles column
db = open (infile_name, "r")

for line in db:
	if line.startswith(">") == False:
		line = line.strip().split("\t")[-1]
		multispecies = False
		if "MULTISPECIES:" in line:
			multispecies = True
			line = line.replace("MULTISPECIES: ","").strip()

		db_line_counter += 1
		splitline = line.split("[",1)
		# ID, the hit returned in DIAMOND results
		db_id = str(splitline[0].split()[0]) # [1:] # Not needed anymore

		# name and functional description
		db_entry = line.split("[", 1)
		# 1 ['WP_168619497.1 ATP-binding cassette domain-containing protein, partial ', 'Rhizobium ruizarguesonis]']
		db_entry = db_entry[0].split(" ", 1)
		# 2 ['WP_168619497.1', 'ATP-binding cassette domain-containing protein, partial ']
		db_entry = db_entry[1][:-1]
		# 3 ATP-binding cassette domain-containing protein, partial

		# organism name
		if line.count("[") != 1:
			splitline = line.split("[")

			db_org = splitline[line.count("[")].strip()[:-1]
			if db_org[0].isdigit():
				split_db_org = db_org.split()
				try:
					if split_db_org[1] == "sp.":
						db_org = split_db_org[0] + " " + split_db_org[1] + " " + split_db_org[2]
					else:
						db_org = split_db_org[1] + " " + split_db_org[2]
				except IndexError:
					try:
						db_org = split_db_org[1]
					except IndexError:
						db_org = splitline[line.count("[")-1]
						if db_org[0].isdigit():
							split_db_org = db_org.split()
							db_org = split_db_org[1] + " " + split_db_org[2]
		else:
			db_org = line.split("[", 1)
			db_org = db_org[1].split()
			try:
				db_org = str(db_org[0]) + " " + str(db_org[1])
			except IndexError:
				db_org = line.strip().split("[", 1)
				# Adding sp. to group with other unknown species of the same genus
				db_org = db_org[1][:-1] + " sp."
				# Multispecies not counted in the error counter
				if not multispecies:
					db_error_counter += 1

		db_org = re.sub('[^a-zA-Z0-9-_*. ]', '', db_org)

		# add to dictionaries
		if "-F" in sys.argv:
			db_func_dictionary[db_id] = db_entry
		if "-O" in sys.argv:
			db_org_dictionary[db_id] = db_org
		if "-SO" in sys.argv:
			if target_org in db_org:
				db_SO_dictionary[db_id] = db_entry
	else:
		print("ELSE SELSLELELSL")

print ("\nSuccess!")
print ("Number of lines: " + str(db_line_counter))
print ("Number of errors: " + str(db_error_counter))

# condensing down the identical matches
condensed_RefSeq_hit_db = {}

for entry in RefSeq_hit_count_db.keys():
	try:
		if "-O" in sys.argv:
			org = db_org_dictionary[entry]
		if "-F" in sys.argv:
			org = db_func_dictionary[entry]
		if org in condensed_RefSeq_hit_db.keys():
			condensed_RefSeq_hit_db[org] += RefSeq_hit_count_db[entry]
		else:
			condensed_RefSeq_hit_db[org] = RefSeq_hit_count_db[entry]
	except KeyError:
		print ("KeyError:\t" + entry)
		continue

if "SO" in sys.argv:
	condensed_RefSeq_SO_hit_db = {}

	for entry in RefSeq_hit_count_db.keys():
		if entry in db_SO_dictionary.values():
			org = db_SO_dictionary[entry]
			if org in condensed_RefSeq_SO_hit_db.keys():
				condensed_RefSeq_SO_hit_db[org] += RefSeq_hit_count_db[entry]
			else:
				condensed_RefSeq_SO_hit_db[org] = RefSeq_hit_count_db[entry]


# dictionary output and summary
print ("\nDictionary database assembled.")
print ("Number of errors: " + str(db_error_counter))

if "-O" in sys.argv:
	print ("\nTop ten organism matches:")
if "-F" in sys.argv:
	print ("\nTop ten function matches:")
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda kv: -kv[1])[:10]:
	try:
		print (str(v) + "\t" + k )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue

# creating the outfiles
if "-O" in sys.argv:
	outfile_name = infile_name[:-4] + "_organism.tsv"
if "-F" in sys.argv:
	outfile_name = infile_name[:-4] + "_function.tsv"
if "=SO" in sys.argv:
	target_org_outfile = open(infile_name[:-4] + "_" + target_org + ".tsv", "w")

outfile = open (outfile_name, "w")

# writing the output
error_counter = 0
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda kv: -kv[1]):
	try:
		q = v * 100 / float(line_counter)
		outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
	except KeyError:
		outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue

# writing the output if optional specific organism flag is active
if "-SO" in sys.argv:
	for k, v in sorted(condensed_RefSeq_SO_hit_db.items(), key=lambda kv: -kv[1]):
		try:
			q = v * 100 / float(line_counter)
			target_org_outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
		except KeyError:
			target_org_outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
			error_counter += 1
			continue

print ("\nAnnotations saved to file: '" + outfile_name + "'.")
print ("Number of errors: " + str(error_counter))

db.close()
outfile.close()
