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
# -D		database		specifies a reference database to search against
#								for results
# -O		organism		returns organism results
# -F		function		returns functional results
# -R		reference		returns reference IDs in results
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

t0 = time.time()

# checking for an option (organism or function) to be specified
if "-O" not in sys.argv:
	if "-F" not in sys.argv:
		if "-R" not in sys.argv:
			sys.exit("WARNING: need to specify either organism results (with -O flag in command), reference IDs (with -R flag in command), or functional results (with -F flag in command).")

# loading starting file
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

# checking to make sure database is specified
if "-D" in sys.argv:
	db_name = string_find("-D")
else:
	sys.exit( "No database file indicated; skipping database search step.")


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
	if line_counter % 1000000 == 0:
		t99 = time.time()
		print (str(line_counter)[:-6] + "M lines processed so far in " + str(t99-t0) + " seconds.")

	unique_seq_db[splitline[0]] = 1

	try:
		RefSeq_hit_count_db[splitline[1]] += 1
	except KeyError:
		RefSeq_hit_count_db[splitline[1]] = 1
		continue

t1 = time.time()

print ("\nAnalysis of " + infile_name + " complete.")
print ("Number of total lines: " + str(line_counter))
print ("Number of unique sequences: " + str(len(unique_seq_db)))
print ("Time elapsed: " + str(t1-t0) + " seconds.")

infile.close()

# time to search for these in the reference database
db = open (db_name, "r")

print ("\nStarting database analysis now.")

t2 = time.time()

# optional outfile of specific organism results
if "-SO" in sys.argv:
	target_org = string_find("-SO")
	db_SO_dictionary = {}

# building a dictionary of the reference database
if "-F" in sys.argv:
	db_func_dictionary = {}
if "-O" in sys.argv:
	db_org_dictionary = {}
if "-R" in sys.argv:
	db_ref_dictionary = {}
db_line_counter = 0
db_error_counter = 0

for line in db:
	if line.startswith(">") == True:
		db_line_counter += 1
		splitline = line.split()

		# ID, the hit returned in DIAMOND results
		db_id = str(splitline[0])[1:]

		# name and functional description
		db_entry = line.split("[", 1)
		db_entry = db_entry[0].split(" ", 1)
		db_entry = db_entry[1][:-1]

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
				db_org = db_org[1][:-1]
				db_error_counter += 1

		db_org = re.sub('[^a-zA-Z0-9-_*. ]', '', db_org)

		# add to dictionaries
		if "-F" in sys.argv:
			db_func_dictionary[db_id] = db_entry
		if "-O" in sys.argv:
			db_org_dictionary[db_id] = db_org
		if "-R" in sys.argv:
			db_ref_dictionary[db_id] = db_id
		if "-SO" in sys.argv:
			if target_org in db_org:
				db_SO_dictionary[db_id] = db_entry

		# line counter to show progress
		if db_line_counter % 1000000 == 0:							# each million
			t95 = time.time()
			print (str(db_line_counter)[:-6] + "M lines processed so far in " + str(t95-t2) + " seconds.")

t3 = time.time()

print ("\nSuccess!")
print ("Time elapsed: " + str(t3-t2) + " seconds.")
print ("Number of lines: " + str(db_line_counter))
print ("Number of database entries with only genus name (no species): " + str(db_error_counter))

# condensing down the identical matches
condensed_RefSeq_hit_db = {}

for entry in RefSeq_hit_count_db.keys():
	try:
		if "-O" in sys.argv:
			org = db_org_dictionary[entry]
		if "-F" in sys.argv:
			org = db_func_dictionary[entry]
		if "-R" in sys.argv:
			org = org + "\t" + entry
		if org in condensed_RefSeq_hit_db.keys():
			condensed_RefSeq_hit_db[org] += RefSeq_hit_count_db[entry]
		else:
			condensed_RefSeq_hit_db[org] = RefSeq_hit_count_db[entry]
	except KeyError:
		print ("KeyError:\t" + entry)
		continue

if "-SO" in sys.argv:
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
print ("Time elapsed: " + str(t3-t2) + " seconds.")
# print ("Number of errors - IDs not found in reference database: " + str(db_error_counter))

if "-O" in sys.argv:
	print ("\nTop ten organism matches:")
if "-F" in sys.argv:
	print ("\nTop ten function matches:")
if "-R" in sys.argv:
	print ("\nTop ten reference IDs:")
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
if "-R" in sys.argv:
	outfile_name = infile_name[:-4] + "_referenceIDs.tsv"
if "-SO" in sys.argv:
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

db.close()
outfile.close()

print ("\nAnnotations saved to file: '" + outfile_name + "'.")
print ("Number of errors in writing annotations to outfile: " + str(error_counter))
