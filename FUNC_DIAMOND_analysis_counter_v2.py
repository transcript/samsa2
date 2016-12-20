#!/usr/bin/env Python

# DIAMOND_analysis_counter.py
# Created 8/16/2016, by Sam Westreich, swestreich@gmail.com, github.com/transcript
# Last modified 11/09/2016

# Purpose: This parses through the results file from a DIAMOND annotation run (in BLAST
# m8 format) to get the results into something a bit more compressed and readable.

# imports
import operator, sys, time, gzip, re

# usage statement
# TO COME LATER ***************************
#	-i		input file					REQUIRED
#	-d		database file				REQUIRED
#	-o		organism-specific mode		OPTIONAL

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem

t0 = time.clock()

# loading starting file
if "-i" in sys.argv:
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
temp_outfile_name = infile_name[:-5] + ".out"
temp_outfile  = open (temp_outfile_name, "w")
for k, v in sorted(RefSeq_hit_count_db.items(), key=lambda (k,v): -v):
	q = v * 100 / float(line_counter)
	temp_outfile.write (str(q) + "%\t" + str(v) + "\t" + k + "\n")

print "Results are saved in " + temp_outfile_name + " for quick resuming."

temp_outfile.close()

print ("\nTop ten matches:")
for k, v in sorted(RefSeq_hit_count_db.items(), key=lambda (k,v): -v)[:10]:
	q = v * 100 / float(line_counter)
	print (str(q) + "%\t" + str(v) + "\t" + k )
	
# time to search for these in the reference database
if "-d" in sys.argv:
	db_name = string_find("-d")
else:
	sys.exit( "No database file indicated; skipping database search step.")

db = open (db_name, "r")

print "\nStarting to build database dictionary now."

t2 = time.clock()

# building a dictionary of the reference database
db_func_dictionary = {}
db_org_dictionary = {}
db_line_counter = 0
db_error_counter = 0

# organism-specific functional results
org_specific_flag = False
if "-o" in sys.argv:
	organism = string_find("-o")
	org_specific_flag = True

# searching 
for line in db:
	if line.startswith(">") == True:	
		db_line_counter += 1
		splitline = line.split("  ")
		
		#organism-specific filtering is here
		if org_specific_flag == True:
			if organism not in line:
				continue
		
		# ID, the hit returned in DIAMOND results
		db_id = str(splitline[0])[1:]
		
		# name and functional description
		db_entry = re.split("\[\d\d", line)
		db_entry1 = db_entry[0].split(" ", 1)
		db_entry2 = db_entry1[1][1:-1]
		
		# function name
		try:
			db_func = db_entry2.split(" ", 1)[1]
			try:
				if db_func[0] == " ":
					db_func = db_func[1:]
			except IndexError:
				db_entry = line.split("  ")
				db_func = db_entry[2]
				continue
		except IndexError:
			print "FUNCTION ERROR ON:\t" + line
			db_func = line.split("  ")[2]
			continue
		
#		# add to dictionaries		
		db_func_dictionary[db_id] = db_func
				
		if db_line_counter % 1000000 == 0:									# each million
			t99 = time.clock()
			print str(db_line_counter) + " lines processed so far in " + str(t99-t2) + " seconds."
		
t3 = time.clock()

print "\nSuccess!"
print "Time elapsed: " + str(t3-t2) + " seconds."
print "Number of lines: " + str(db_line_counter)
print "Number of errors: " + str(db_error_counter)

# condensing down the identical matches
print "\nStarting to form the output file now."
condensed_RefSeq_hit_db = {}

condensed_RefSeq_func_hit_db = {}
for entry in RefSeq_hit_count_db.keys():
	try:
		func = db_func_dictionary[entry]
	#	try:
	#		simplified_func = func.split(" ", 1)[1]
	#		if simplified_func.strip() == "":
	#			simplified_func = func
	#	except IndexError:
	#		print func
	#		simplified_func = func
	#		continue
		simplified_func = func
		if simplified_func in condensed_RefSeq_func_hit_db.keys():
			condensed_RefSeq_func_hit_db[simplified_func] += RefSeq_hit_count_db[entry]
		else:
			condensed_RefSeq_func_hit_db[simplified_func] = RefSeq_hit_count_db[entry]
	except KeyError:
		continue

# dictionary output and summary
print "\nDictionary database assembled."
print "Time elapsed: " + str(t3-t2) + " seconds."
print "Number of errors: " + str(db_error_counter)

print "\nTop ten function matches:"
for k, v in sorted(condensed_RefSeq_func_hit_db.items(), key=lambda (k,v): -v)[:10]:
	try:
		print (str(v) + "\t" + k ) #db_func_dictionary[k] )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue
		
# creating the outfiles
if "-o" in sys.argv:
	func_outfile_name = infile_name[:-5] + "_" + organism + "_function.tsv"
else:
	func_outfile_name = infile_name[:-5] + "_function.tsv"

func_outfile = open (func_outfile_name, "w")

# output for functional annotations
error_counter = 0
for k, v in sorted(condensed_RefSeq_func_hit_db.items(), key=lambda (k,v): -v):
	try:
		q = v * 100 / float(line_counter)
		func_outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n") #db_func_dictionary[k] + "\n")
	except KeyError:
		func_outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue

print "\nFunctional annotations saved to file: '" + func_outfile_name + "'."
print "Number of errors: " + str(error_counter)

db.close()
func_outfile.close()