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
# raw_read_counter.py
# Created 6/19/2017, this version edited 6/19/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# PURPOSE: When running DESeq stats on files, the total number of input
# sequences needs to be included.  This program generates an outfile with
# the number of reads per sequence file, for later import into the SAMSA
# pipeline.
##########################################################################
#
# USAGE:
#
# python raw_read_counter.py -I $infile -O $stats_outfile
#
##########################################################################

# imports
import sys, time, gzip

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		this_elem = elem
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem

if "-I" not in sys.argv:
	sys.exit("WARNING: Specify the infile with the '-I' flag.")
else:
	infile_name = string_find("-I")
if "-O" not in sys.argv:
	sys.exit("WARNING: Specify the outfile with the '-O' flag.")
else:
	outfile_name = string_find("-O")

# setting variables
read_counter = 0

# reading through infile
infile = open(infile_name, "r")

for line in infile:
	read_counter += 1

infile.close()

read_count = read_counter / 4 # to account for fastq files

# writing to outfile
outfile = open(outfile_name, "a")

file_name = infile_name.split(".cleaned.fastq")[0]

if "/" in file_name:
	file_name_cut = file_name.split("/")[-1]
	file_name = file_name_cut

outfile.write(file_name + "\t" + str(read_count) + "\n")

outfile.close()
