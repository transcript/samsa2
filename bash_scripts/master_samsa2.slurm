#!/bin/bash

#SBATCH --time=0-23:59
#SBATCH --mem=100000
#SBATCH --account=def-bruneaua
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48

module load r/3.5.0


# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM
####################################################################
#set -x #echo
#
# master_script_for_sample_files.bash
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified by Michelle Treiber on August 9, 2017
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#	0.1	Create the checkpoint
#	0.2	Gunzip the raw sequences
#	1.	Merging with PEAR, if applicable
#	2.	Read cleaning with Trimmomatic
#	3.	rRNA removal with SortMeRNA
#	4.	Annotation using DIAMOND against Subsys database
#	5.	Aggregation using analysis_counter.py (Subsys)
#   	4.1  	Annotation using DIAMOND against Refseq database
#   	5.1	Aggregation using analysis_counter.py (Refseq)
#	6.	Running R scripts to get DESeq statistical analysis. (NOT RUN)
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set 
# up SAMSA2 dependencies and download full databases.
#
####################################################################
#
echo -e "NOTE: Before running this script, user must run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - Set pathway for starting_location to location of samsa2 GitHub download
#
# 	0. Set starting location:
starting_location=/home/renaut/scratch/ATRAPP_Champlain_2016_Metat/metagenome_St2

#	0.5 Sequences location (note that sequences actually need to be stored in directory called $sequence_location/hiseq_raw
sequence_location=/home/renaut/scratch/ATRAPP_Champlain_2016_Metat/sequences/metagenome/St2

#	1. PEAR
pear_location=/home/renaut/bin

# 	2. Trimmomatic
trimmomatic_location=/home/renaut/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic-0.36

#	3. SortMeRNA
sortmerna_location=/home/renaut/bin

#	4. DIAMOND
diamond_database="/home/renaut/scratch/blast_database/refseq/RefSeq_bac"
diamond_subsys_db="/home/renaut/scratch/blast_database/subsys/subsys_db"
diamond_location=/home/renaut/bin

#	5. Aggregation
python_programs=/home/renaut/samsa2.master/samsa2/python_scripts
RefSeq_db="/home/renaut/scratch/blast_database/refseq/RefSeq_bac.fa"
Subsys_db="/home/renaut/scratch/blast_database/subsys/subsys_db.fa"

#	6. R scripts and paths
export R_LIBS="/home/renaut/samsa2.master/samsa2/R_scripts/packages"
R_programs=/home/renaut/samsa2.master/samsa2/R_scripts

#	threads
threads="48"

##################
#STEP 0.1: create/read checkpoint
cd $starting_location

printf "\nStep 0.1: Checking for the presence of the checkpoint file\n"
if [ ! -f "pipeline/checkpoints" ]
        then
                printf "\tThe file checkpoints does not exist, we will create it.\n"
        touch "pipeline/checkpoints"
else
        printf "\tThe file checkpoints already exist\n"
fi

####################################################################
#
#STEP 0.2: GUNZIP PROCESS. 
#
#
Step=$(grep "GUNZIP" pipeline/checkpoints)
if [ "${Step}" != "GUNZIP" ]
	then
		for file in $sequence_location/hiseq_raw/*gz
			do
				file_out=`echo $file | awk -F ".gz" '{print $1}'`	
        			gunzip $file -cd >$file_out

				if [ $? -ne 0 ]
                        		then
                                	printf "\t!!! There is a problem in the gunzip step"
                                	exit 1
                		fi
			done
		
		mkdir $sequence_location/hiseq_unzip
		mv $sequence_location/hiseq_raw/*fastq $sequence_location/hiseq_unzip/.
		printf "GUNZIP\n" >>pipeline/checkpoints
                echo "STEP 0.1 DONE:GUNZIP OF PAIRED-END FILES USING GUNZIP"
                date
else
        printf  "\tThe variable GUNZIP is in the checkpoint file. STEP 0.1 will then be passed\n"
fi

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#		Example: control_1.R1.fastq
#				 control_1.R2.fastq
#
# Note: if using single-end sequencing, skip this step (comment out).
# Note: if performing R analysis (step 6), be sure to name files with 
# 	the appropriate prefix ("control_$file" and "experimental_$file")!



Step=$(grep "MERGING" pipeline/checkpoints)
if [ "${Step}" != "MERGING" ]
	then
		mkdir $starting_location/step_1_merging
		touch $starting_location/step_1_merging/pear_log
		for file in $sequence_location/hiseq_unzip/*_R1.fastq
		do
			file1=$file
			file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'`
			out_path=`echo $file | awk -F "R1" '{print $1 "merged"}'`
			out_name=`echo ${out_path##*/}`
			$pear_location/pear -f $file1 -j $threads -r $file2 -o $out_name 1>>$starting_location/step_1_merging/pear_log
		
			if [ $? -ne 0 ]
                        	then
                                	printf "\t!!! There is a problem in the merging step"
                                	exit 1
                	fi
		done
		
		###I will concatenate forward and reverse complement with 20Ns in the middle.
		$R_programs/combining_umerged.R $starting_location
		
		printf "MERGING\n" >>pipeline/checkpoints
		mkdir $starting_location/step_1_merging/
		mv $starting_location/*assembled2.fastq $starting_location/step_1_merging/
		
		rm $starting_location/*fastq
		rm $sequence_location/hiseq_unzip/*fastq

		echo "STEP 1 DONE:MERGING OF PAIRED-END FILES USING PEAR"
		date        
else
        printf  "\tThe variable MERGING is in the checkpoint file. STEP 1 will then be passed\n"
fi

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $starting_location/step_1_output/ folder!
Step=$(grep "TRIMMING" pipeline/checkpoints)
if [ "${Step}" != "TRIMMING" ]
        then
		mkdir $starting_location/step_2_trimming
		touch $starting_location/step_2_trimming/trimmomatic_log

	for file in $starting_location/step_1_merging/*assembled2.fastq
		do
		shortname=`echo $file | awk -F "assembled2.fastq" '{print $1 "cleaned.fastq"}'`

		java -jar $trimmomatic_location/trimmomatic-0.36.jar SE -threads $threads -phred33 $file $shortname ILLUMINACLIP:$trimmomatic_location/adapters/all_primers.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:59 2>>$starting_location/step_2_trimming/trimmomatic_log
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the trimmmomatic step"
				exit 1
		fi
		
		done

	mkdir $starting_location/step_2_trimming/
	mv $starting_location/step_1_merging/*cleaned.fastq $starting_location/step_2_trimming/
	
	rm step_1_merging/*assembled2.fastq
	
	echo "STEP 2 DONE:CLEANING FILES WITH TRIMMOMATIC"
	date
	
	printf "TRIMMING\n" >>pipeline/checkpoints
else
        printf  "\tThe variable TRIMMING is in the checkpoint file. STEP 2 will then be passed\n"
fi

####################################################################
#
# STEP 2.5: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.
Step=$(grep "RAW" pipeline/checkpoints)
if [ "${Step}" != "RAW" ]
        then
		if [ -f $starting_location/step_2_trimming/raw_counts.txt ]
			then
				rm $starting_location/step_2_trimming/raw_counts.txt
				touch $starting_location/step_2_trimming/raw_counts.txt
		else
			touch $starting_location/step_2_trimming/raw_counts.txt
		fi	

		for file in $starting_location/step_2_trimming/*cleaned.fastq
		do
			python $python_programs/raw_read_counter.py -I $file -O $starting_location/step_2_trimming/raw_counts.txt
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the raw sequences step"
				exit 1
		fi
		done
		echo "STEP 2.5:GETTING RAW SEQUENCES COUNTS" 
		date

        printf "RAW\n" >>pipeline/checkpoints
else
        printf  "\tThe variable RAW is in the checkpoint file. STEP 2.9 will then be passed\n"
fi

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

Step=$(grep "RIBODEPLETION" pipeline/checkpoints)
if [ "${Step}" != "RIBODEPLETION" ]
        then
        mkdir $starting_location/step_3_ribodepletion
		touch $starting_location/step_3_ribodepletion/sortmerna_log

		for file in $starting_location/step_2_trimming/*cleaned.fastq
		do
			shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`
			$sortmerna_location/sortmerna -a $threads --ref $sortmerna_location/../sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/../sortmerna-2.1/index/silva-bac-16s --reads $file --aligned $file.ribosomes --other $shortname --fastx --log -v 1>>$starting_location/step_3_ribodepletion/sortmerna_log
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the sortmerna step"
				exit 1
		fi
		done

		mkdir $starting_location/step_3_ribodepletion/
		mv $starting_location/step_2_trimming/*ribodepleted* $starting_location/step_3_ribodepletion/
		
		rm step_2_trimming/*cleaned.fastq
		rm step_2_trimming/*ribosomes.fastq
		
		echo "STEP 3 DONE:REMOVING RIBOSOMAL READS WITH SORTMERNA"
		date
	
        	printf "RIBODEPLETION\n" >>pipeline/checkpoints
else
        	printf  "\tThe variable RIBODEPLETION is in the checkpoint file. STEP 3 will then be passed\n"
fi

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS
date
Step=$(grep "ANNOTATION_SUBSYS" pipeline/checkpoints)
if [ "${Step}" != "ANNOTATION_SUBSYS" ]
	then
		mkdir $starting_location/step_4_annotation
		touch $starting_location/step_4_annotation/diamond_log
		for file in $starting_location/step_3_ribodepletion/*ribodepleted.fastq
		do
			shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
			#removed the diamond sensitive option
			$diamond_location/diamond blastx --db $diamond_subsys_db -q $file -a $file.Subsys -p $threads -t ./ -k 1 1>>$starting_location/step_4_annotation/diamond_log
			$diamond_location/diamond view --daa $file.Subsys.daa -o $shortname -p $threads -f tab 1>>$starting_location/step_4_annotation/diamond_log
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the diamond subsys step"
				exit 1
		fi
		done
		
		mkdir $starting_location/step_4_annotation/daa_binary_files/
	
		mv $starting_location/step_3_ribodepletion/*subsys_annotated* $starting_location/step_4_annotation/
		mv $starting_location/step_3_ribodepletion/*.daa $starting_location/step_4_annotation/daa_binary_files/

		echo "STEP 4 DONE: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS"
        	date

		printf "ANNOTATION_SUBSYS\n" >>pipeline/checkpoints
else
       	 printf  "\tThe variable ANNOTATION_SUBSYS is in the checkpoint file. STEP 4 will then be passed\n"
fi
##################################################################
#
# STEP 5: PYTHON SUBSYSTEMS ANALYSIS COUNTER
Step=$(grep "AGGREGATION_SUBSYS" pipeline/checkpoints)
if [ "${Step}" != "AGGREGATION_SUBSYS" ]
	then
		mkdir $starting_location/step_5_aggregation/
		for file in $starting_location/step_4_annotation/*subsys_annotated
		do
			python $python_programs/DIAMOND_subsystems_analysis_counter.py -I $file -D $Subsys_db -O $file.hierarchy -P $file.receipt
	
			# This quick program reduces down identical hierarchy annotations
			python $python_programs/subsys_reducer.py -I $file.hierarchy
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the aggregation subsys step"
				exit 1
		fi
		done

		mkdir $starting_location/step_5_aggregation/Subsystems_results/
		mkdir $starting_location/step_5_aggregation/Subsystems_results/receipts/
		mv $starting_location/step_4_annotation/*.reduced $starting_location/step_5_aggregation/Subsystems_results/
		mv $starting_location/step_4_annotation/*.receipt $starting_location/step_5_aggregation/Subsystems_results/receipts/
		rm $starting_location/step_4_annotation/*.hierarchy

		echo "STEP 5 DONE:PYTHON SUBSYSTEMS ANALYSIS COUNTER"
		date

        	printf "AGGREGATION_SUBSYS\n" >>pipeline/checkpoints
else
	        printf  "\tThe variable AGGREGATION_SUBSYS is in the checkpoint file. STEP 5. will then be passed\n"
fi


####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.
Step=$(grep "ANNOTATION_REFSEQ" pipeline/checkpoints)
if [ "${Step}" != "ANNOTATION_REFSEQ" ]
        then
		touch $starting_location/step_4_annotation/diamond_log

		for file in $starting_location/step_3_ribodepletion/*ribodepleted.fastq
			do
				shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
				#removed the diamond sensitive option
				$diamond_location/diamond blastx --db $diamond_database -q $file -p $threads -a $file.RefSeq -t ./ -k 1 1>>$starting_location/step_4_annotation/diamond_log
				$diamond_location/diamond view --daa $file.RefSeq.daa -o $shortname -p $threads -f tab 1>>$starting_location/step_4_annotation/diamond_log
			if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the diamond refseq step"
				exit 1
			fi
			done

		mv $starting_location/step_3_ribodepletion/*annotated* $starting_location/step_4_annotation/
		mv $starting_location/step_3_ribodepletion/*.daa $starting_location/step_4_annotation/daa_binary_files/

		echo "STEP 4.1 DONE:ANNOTATING WITH DIAMOND AGAINST REFSEQ"
		date

        	printf "ANNOTATION_REFSEQ\n" >>pipeline/checkpoints
else
	        printf  "\tThe variable ANNOTATION_REFSEQ is in the checkpoint file. STEP 4.1 will then be passed\n"
fi

####################################################################
#
# STEP 5.1: AGGREGATING WITH ANALYSIS_COUNTER
Step=$(grep "AGGREGATION_REFSEQ" pipeline/checkpoints)
if [ "${Step}" != "AGGREGATION_REFSEQ" ]
	then
		for file in $starting_location/step_4_annotation/*RefSeq_annotated*
		do
			python $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
			python $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem in the aggregation refseq step"
					exit 1
			fi
		done

		mkdir $starting_location/step_5_aggregation/RefSeq_results/
		mkdir $starting_location/step_5_aggregation/RefSeq_results/org_results/
		mkdir $starting_location/step_5_aggregation/RefSeq_results/func_results/
		mv $starting_location/step_4_annotation/*organism.tsv $starting_location/step_5_aggregation/RefSeq_results/org_results/
		mv $starting_location/step_4_annotation/*function.tsv $starting_location/step_5_aggregation/RefSeq_results/func_results/

		echo "STEP 5.1 DONE: AGGREGATING WITH ANALYSIS_COUNTER"
		date 

        	printf "AGGREGATION_REFSEQ\n" >>pipeline/checkpoints
else
	        printf  "\tThe variable AGGREGATION_REFSEQ is in the checkpoint file. STEP 5.1 will then be passed\n"
fi



##################################################################
#
# At this point, all the results files are ready for analysis using R.
# This next step performs basic DESeq2 analysis of the RefSeq organism, function,
# and Subsystems annotations.
#
# More complex R analyses may be performed using specific .sh analysis scripts.
#
# STEP 6: R ANALYSIS
# Note: For R to properly identify files to compare/contrast, they must include 
# the appropriate prefix (either "control_$file" or experimental_$file")!

#/home/apps/Logiciels/R/3.4.3-gcc/bin/Rscript $R_programs/run_DESeq_stats.R -I $starting_location/step_5_aggregation/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab -R $starting_location/step_2_trimming/raw_counts.txt
#/home/apps/Logiciels/R/3.4.3-gcc/bin/Rscript $R_programs/run_DESeq_stats.R -I $starting_location/step_5_aggregation/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab -R $starting_location/step_2_trimming/raw_counts.txt
#/home/apps/Logiciels/R/3.4.3-gcc/bin/Rscript $R_programs/Subsystems_DESeq_stats.R -I $starting_location/step_5_aggregation/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1 -R step_2_trimming/raw_counts.txt

echo "STEP 6 DONE: R ANALYSIS"
echo "Master bash script finished running at: "; date
exit
####################################################################

