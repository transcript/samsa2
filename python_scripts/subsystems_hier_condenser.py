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
# subsystems_hier_condenser.py
# Created 2/01/2017, this version created 2/01/2017
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
# -H		hierarcy level	level of output requested (3, 2, or 1)
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

# loading starting file, hierarchy level
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

if "-H" in sys.argv:
	hier_lvl = string_find("-H")
	if int(hier_lvl) not in [1,2,3]:
		sys.exit("WARNING: requested hierarchy level must be 1, 2, or 3.")
else:
	sys.exit("WARNING: Specify requested hierarchy level using '-H' flag.")

infile = open (infile_name, "r")

hier_count_dic = {}
hier_percent_dic = {}
line_counter = 0

for line in infile:
	line_counter += 1
	splitline = line.strip().split("\t")
	if "NO HIERARCHY" in line:
		continue

	# hierarchy level 1 - most broad	
	elif int(hier_lvl) == 1:
		try:
			hier_count_dic[splitline[5]] += int(splitline[1])
			hier_percent_dic[splitline[5]] += float(splitline[0])
		except KeyError:
			hier_count_dic[splitline[5]] = int(splitline[1])
			hier_percent_dic[splitline[5]] = float(splitline[0])

	# hierarchy level 2 - middle level
	elif int(hier_lvl) == 2:
		try:
			hier_count_dic[splitline[4]] += int(splitline[1])
			hier_percent_dic[splitline[4]] += float(splitline[0])
		except KeyError:
			hier_count_dic[splitline[4]] = int(splitline[1])
			hier_percent_dic[splitline[4]] = float(splitline[0])

	# hierarchy level 3 - most specific, besides function
	elif int(hier_lvl) == 3:
		try:
			hier_count_dic[splitline[3]] += int(splitline[1])
			hier_percent_dic[splitline[3]] += float(splitline[0])
		except KeyError:
			hier_count_dic[splitline[3]] = int(splitline[1])
			hier_percent_dic[splitline[3]] = float(splitline[0])

infile.close()

# outfile time
outfile = open(infile_name + ".level_" + hier_lvl, "w")
for key, value in sorted(hier_count_dic.items(), key=lambda(key, value): -value):
	outfile.write(str(hier_percent_dic[key]) + "\t" + str(value) + "\t" + key + "\n")
	
print "Converted data saved in " + infile_name + ".level_" + hier_lvl
outfile.close()

print "\nTop ten hierarchy matches at new level:"
for k, v in sorted(hier_count_dic.items(), key=lambda (k,v): -v)[:10]:
	try:
		print (str(v) + "\t" + k )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue
	
