#!/bin/bash

####################################################################
#
# DIAMOND_example_script.bash
# Created December 2016 by Sam Westreich, github.com/transcript
#
####################################################################
#
# This bash script shows how to structure commands for the following
# purposes:
#
# 	1. Converting a reference database to DIAMOND-searchable format;
#	2. Performing an annotation search, running an input file against
#	   a DIAMOND-searchable database;
#	3. Converting the results of an annotation search to BLAST m8
#      format.
#
####################################################################

# setting variables:

database = "starting_reference_database_name.fa"
diamond_database = "starting_reference_database_name.daa"
filename = "starting_fastq_input_file.fastq"
diamond_location = /path/to/diamond/program/
diamond_output = "name_of_DIAMOND_search_output.daa"
final_output = "name_of_final_BLAST_m8_results.m8"

####################################################################

# DATABASE CREATION EXAMPLE

echo "Creating database from " $database
date

$diamond_location/diamond makedb --in $database --db $database 

# explanation of settings:
#
#	--in	starting file that will be used to create DIAMOND-searchable db
#	--db	name of created DIAMOND database (automatically appended with .daa)

####################################################################

# ANNOTATION SEARCH EXAMPLE

echo "Performing annotation search on " $filename " against " $database
date

$diamond_location/diamond blastx --db $diamond_database -q $filename -a $diamond_output -t ./ -k 1

# explanation of settings:
#
#	--db	sets database (must be in DIAMOND-readable form)
#	-q		query file name
#	-a		name of results file (in DIAMOND format)
#	-t		sets temporary directory locations
#	-k		number of hits above cutoff to return (if run without specifying,
#			default is 25)

####################################################################

# CONVERTING RESULTS EXAMPLE

echo "Converting file " $diamond_output " to readable format"
date

$diamond_location/diamond view --daa $diamond_output -o $final_output -f tab

# explanation of settings:
#
#	--daa	name of output file in .daa format, from the annotation search
#	-o		output file name
#	-f		separator of values in final output

####################################################################