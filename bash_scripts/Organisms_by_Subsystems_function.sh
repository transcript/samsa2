#!/bin/bash
#SBATCH --mem=32000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user=stwestreich@ucdavis.edu
#SBATCH --mail-type=END

####################################################################
#
# Organisms_by_Subsystems_function.sh
# Created June 17, 2017 by Sam Westreich, github.com/transcript
# This version modified June 17, 2017
#
####################################################################
#
# Retrieving all organism reads from RefSeq that match a Subsystems 
# functional annotation.
#
####################################################################
#
# VARIABLES

python_programs=/share/milklab/sam/python_scripts

specific_function="Carbohydrates"

RefSeq_db="/share/milklab/sam/databases/bct"

RefSeq_results_location=/share/milklab/sam/test_files/output_diamond
Subsys_results_location=/share/milklab/sam/test_files/Subsys_results/receipts
base_directory=/share/milklab/sam/test_files

export R_LIBS="/share/milklab/sam/R_scripts/packages"
R_programs=/share/milklab/sam/R_scripts

#
####################################################################
#
# PREREQUISITES
#
# For this program to work, you should have already run:
#   * DIAMOND_analysis_counter.py against the RefSeq results
#   * DIAMOND_Subsystems_analysis_counter.py against the Subsystems
#		results (with -P flag to obtain a .receipt file)
#
####################################################################
#
# STEP 1 - Python, mapping the Subsystems results matching the 
#	"specific function" (in variables above) to hits to get the read IDs.

for file in $Subsys_results_location/*.receipt
do
	# gets file name without directory name
	file_name=`echo $file | awk -F "/" '{print $NF}'`
	
	shortname=`echo $file_name | awk -F ".receipt" '{print $1}'`
	
	RefSeq_results_name=$shortname.RefSeq
	output_name=$shortname.$specific_function.org
	
	python $python_programs/Subsys_to_RefSeq_mapper.py -S $file -I $RefSeq_results_location/$RefSeq_results_name -T $specific_function -O $Subsys_results_location/$output_name

done

#
####################################################################
#
# STEP 2 - Python, aggregate these results

for file in $Subsys_results_location/*$specific_function.org
do
	python $python_programs/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
done

#
####################################################################
#
# STEP 2.9 - Moving everything to its own directory

new_dir_suffix="_organisms"
new_dir_name=$specific_function$new_dir_suffix

mkdir $base_directory/$new_dir_name
mkdir $base_directory/$new_dir_name/aggregated

mv $Subsys_results_location/*$specific_function* $base_directory/$new_dir_name
mv $base_directory/$new_dir_name/*organism.tsv $base_directory/$new_dir_name/aggregated/

#
####################################################################
#
# STEP 3 - Run R to get DESeq analysis of the organism results

Rscript $R_programs/run_DESeq_stats.R -I $base_directory/$new_dir_name/aggregated/ -O RefSeq.$specific_function.org_DESeq_results.tab

#
####################################################################
#
# STEP 3.9 - Finishing up, printing exit info

echo "Files have been analyzed, pulling out all organism results to $specific_function hits in the Subsystems results."
echo "The raw results have been saved in $base_directory/$new_dir_name/ ."
echo "The aggregated results have been saved in $base_directory/$new_dir_name/aggregated ."
echo "The DESeq analysis has been saved in the above directory as RefSeq.$specific_function.org_DESeq_results.tab ."


