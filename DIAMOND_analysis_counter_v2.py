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

# checking for an option (organism or function) to be specified
if "-O" not in sys.argv.upper():
	if "-F" not in sys.argv.upper():
		sys.exit("WARNING: need to specify either organism results (with -O flag in command) or functional results (with -F flag in command).")

# loading starting file
if "-i" in sys.argv.lower():
	infile_name = string_find("-i")
else:
	sys.exit ("WARNING: infile must be specified using '-i' flag.")

infile = open (infile_name, "r")

# setting up databases
RefSeq_hit_count_db = {}
unique_seq_db = {}
line_counter = 0

for line in infile:
	line_counter += 1
	splitline = line.split("\t")
	if line_counter % 1000000 == 0:
		t99 = time.clock()
		print str(line_counter)[:-6] + "M lines processed so far in " + str(t99-t0) + " seconds."
	
	unique_seq_db[splitline[0]] = 1

	try:
		RefSeq_hit_count_db[splitline[1]] += 1
	except KeyError:
		RefSeq_hit_count_db[splitline[1]] = 1
		continue

t1 = time.clock()

print "\nAnalysis of " + infile_name + " complete."
print "Number of total lines: " + str(line_counter)
print "Number of unique sequences: " + str(len(unique_seq_db))
print "Time elapsed: " + str(t1-t0) + " seconds."

infile.close()

# temporary outfile
#temp_outfile_name = infile_name[:-5] + ".out"
#temp_outfile  = open (temp_outfile_name, "w")
#for k, v in sorted(RefSeq_hit_count_db.items(), key=lambda (k,v): -v):
#	q = v * 100 / float(line_counter)
#	temp_outfile.write (str(q) + "%\t" + str(v) + "\t" + k + "\n")
#
#print "Results are saved in " + temp_outfile_name + " for quick resuming."

#temp_outfile.close()

#print ("\nTop ten matches:")
#for k, v in sorted(RefSeq_hit_count_db.items(), key=lambda (k,v): -v)[:10]:
#	q = v * 100 / float(line_counter)
#	print (str(q) + "%\t" + str(v) + "\t" + k )
	
# time to search for these in the reference database
if "-d" in sys.argv.lower():
	db_name = string_find("-d")
else:
	sys.exit( "No database file indicated; skipping database search step.")

db = open (db_name, "r")

print "\nStarting database analysis now."

t2 = time.clock()

# building a dictionary of the reference database
if "-F" in sys.argv.upper():
	db_func_dictionary = {}
if "-O" in sys.argv.upper():
	db_org_dictionary = {}
db_line_counter = 0
db_error_counter = 0

for line in db:
	if line.startswith(">") == True:
		db_line_counter += 1
		splitline = line.split("  ")
		
		# ID, the hit returned in DIAMOND results
		db_id = str(splitline[0])[1:]
		
		# name and functional description
		db_entry = line.split("[", 1)
		db_entry = db_entry[0].split(" ", 1)
		db_entry = db_entry[1][1:-1]
		
		# organism name
		if line.count("[") != 1:
			splitline = line.split("[")

			db_org = splitline[line.count("[")].strip()[:-1]
			if db_org[0].isdigit():
				split_db_org = db_org.split()
				try:
					db_org = split_db_org[1] + " " + split_db_org[2]
				except IndexError:
					try:
						db_org = split_db_org[1]
					except IndexError:
						db_org = splitline[line.count("[")-1]
						if db_org[0].isdigit():
							split_db_org = db_org.split()
							db_org = split_db_org[1] + " " + split_db_org[2]
						print line
						print db_org
		else:	
			db_org = line.split("[", 1)
			db_org = db_org[1].split()
			try:
				db_org = str(db_org[1]) + " " + str(db_org[2])
			except IndexError:
				db_org = line.strip().split("[", 1)
				db_org = db_org[1][:-1]
				db_error_counter += 1
		
		db_org = re.sub('[^a-zA-Z0-9-_*. ]', '', db_org)

		# add to dictionaries		
		if "-F" in sys.argv.upper():
			db_func_dictionary[db_id] = db_entry
		if "-O" in sys.argv.upper():
			db_org_dictionary[db_id] = db_org
		
		# line counter to show progress
		if db_line_counter % 1000000 == 0:							# each million
			t95 = time.clock()
			print str(db_line_counter) + " lines processed so far in " + str(t95-t2) + " seconds."
		
t3 = time.clock()

print "\nSuccess!"
print "Time elapsed: " + str(t3-t2) + " seconds."
print "Number of lines: " + str(db_line_counter)
print "Number of errors: " + str(db_error_counter)
#sys.exit()

# condensing down the identical matches
condensed_RefSeq_hit_db = {}

for entry in RefSeq_hit_count_db.keys():
	if "-O" in sys.argv.upper():
		org = db_org_dictionary[entry]
	if "-F" in sys.argv.upper():
		org = db_func_dictionary[entry]
	if org in condensed_RefSeq_hit_db.keys():
		condensed_RefSeq_hit_db[org] += RefSeq_hit_count_db[entry]
	else:
		condensed_RefSeq_hit_db[org] = RefSeq_hit_count_db[entry]

# dictionary output and summary
print "\nDictionary database assembled."
print "Time elapsed: " + str(t3-t2) + " seconds."
print "Number of errors: " + str(db_error_counter)

if "-O" in sys.argv.upper():
	print "\nTop ten organism matches:"
if "-F" in sys.argv.upper():
	print "\nTop ten function matches:"
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda (k,v): -v)[:10]:
	try:
		print (str(v) + "\t" + k )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue

# creating the outfiles
if "-O" in sys.argv.upper():
	outfile_name = infile_name[:-5] + "_organism.tsv"
if "-F" in sys.argv.upper():
		outfile_name = infile_name[:-5] + "_function.tsv"

outfile = open (outfile_name, "w")

# writing the output
error_counter = 0
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda (k,v): -v):
	try:
		q = v * 100 / float(line_counter)
		outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
	except KeyError:
		outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue

print "\Annotations saved to file: '" + outfile_name + "'."
print "Number of errors: " + str(error_counter)

db.close()
outfile.close()
