#!/usr/bin/env Python
import operator, sys

# RefSeq_output_reducer_v1.py
# Created 3/30, Sam Westreich (swestreich@gmail.com)
# This should go through the RefSeq output file (filename.tab.output) and should take
# each line and simplify it down to genus counts, not species.

# String search function
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem

# quiet mode
if "-Q" in sys.argv:
	quiet = True
else:
	quiet = False

if quiet == False:
	print "For usage options, run with flag '-usage'."

# Usage statement
if "-usage" in sys.argv:
	print "USAGE STATEMENT"
	print "-I\tInput file, required"
	print "-O\tOutput file, if not specified will default to 'input.genus'."
	print "-Q\tQuiet mode, hides all messages"
	sys.exit()

# opening input file
try:
	input_file_name = string_find("-I")
except IndexError:
	input_file_name = sys.exit("WARNING: No input file specified.\nSpecify the input file with '-I' flag.")

input_file = open (input_file_name, "r")

# creating output file
if "-O" in sys.argv:
	output_file_name = string_find("-O")
else:
	output_file_name = input_file_name + ".genus"
output_file = open (output_file_name, "w")

# counters
line_counter = 0
total_entries = 0

db = {}

# reading through input file
for line in input_file:
	line_counter += 1
	splitline = line.split("\t")
	Species_name = splitline[2].strip()
	splitname = Species_name.split()
	familyName = splitname[0]
	if familyName in db.keys():
		db[familyName] += int(splitline[1])
	else:
		db[familyName] = int(splitline[1])

	total_entries += int(splitline[1])

# sorting the results from largest to smallest, writing to output file
for k, v in sorted(db.items(), key=lambda (k,v): -v):
	output_file.write (str(v * 100 / float(total_entries)) + "\t" + str(v) + "\t" + k + "\n")

if quiet == False:
	for k, v in sorted(db.items(), key=lambda (k,v): -v)[:10]:
		print (str(v) + "\t" + k)

if quiet == False:
	print ("\nTotal number of entries:\t" + str(total_entries))

input_file.close()
output_file.close()
