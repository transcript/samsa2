# SAMSA, version 2.0

Version 2 of the SAMSA pipeline - faster!  Lighter!  More options!  Less waiting!  Still in pre-alpha testing and design!  

### New in version 2:
* DIAMOND integration, allowing for SAMSA to be run without ever needing an MG-RAST account.
* Option to annotate against custom databases.
* Better, more polished R scripts that can be executed from the command line.
* PCA plots and other graphical outputs.
* More updates to come!

## Quick start
1. Clone or download the following programs:
	1. SAMSA, version 2.0: https://github.com/transcript/samsa_v2
	2. DIAMOND: https://github.com/bbuchfink/diamond
	3. Trimmomatic, a flexible read cleaner: http://www.usadellab.org/cms/?page=trimmomatic
	4. FLASH, if using paired-end data (recommended): https://ccb.jhu.edu/software/FLASH/
2. Preprocess your reads using "SAMSA\_pre\_annotation\_pipeline.py".  
3. Use DIAMOND to annotate your reads against a database of your choosing (note that database must be local and DIAMOND-indexed).  See "example\_DIAMOND\_annotation\_script.bash" for more details.
4. Use "DIAMOND\_analysis\_counter.py" to create a ranked abundance summary of the DIAMOND results from each metatransciptome file.
5. Import these abundance summaries into R and use "run\_DESeq\_stats.R" to determine the most significantly differing features between either individual metatranscriptomes, or control vs. experimental groups.


## SAMSA: Simple Analysis of Metatranscriptomes by Sequence Annotation
Metatranscriptome, RNA-seq data from multiple members of a microbial community, offers incredibly powerful insights into the workings of a complex ecosystem.  RNA sequences are able to not only identify the individual members of a community down to the strain level, but can also provide information on the activity of these microbes at the time of sample collection - something that cannot be determined through other meta- (metagenome, 16S rRNA sequencing) method.  

However, working with metatranscriptome data often proves challenging, given its high complexity and large size.  SAMSA is one of the first bioinformatics pipelines designed with metatranscriptome data specifically in mind.  It accepts raw sequence data in FASTQ form as its input, and performs four phases:

**Preprocessing:** If the sequencing was paired-end, FLASH is used to merge mate pairs.  Trimmomatic is used for the removal of adaptor contamination and low-quality reads.

**Annotation:** Annotation is completed using [DIAMOND, an accelerated BLAST-like sequence aligner.](https://github.com/bbuchfink/diamond)  (Why DIAMOND?  At a standard rate of 10 annotations per second, a standard BLAST approach would take several months to finish - just for a single file!)

**Aggregation:** DIAMOND returns results on a per-read basis, a bit like a ticker tape or a line item receipt.  In the aggregation step, Python scripts condense these line-by-line results to create summary tables.

**Analysis:** R scripts use DESeq to compute most significantly different features between control vs. experimental samples.  These R scripts generate a tabular output with assigned p-values and log2FoldChange scores for each feature.  These 'features' can be either organisms or specific functions.  R can also create graphs showing visual representation of the metatranscriptome(s).

### Individual programs in SAMSA v.2.0 and their functions
For more information, please consult the manual, which goes into more detail on each step in the SAMSA pipeline.

**Preprocessing:**

* SAMSA_pre\_annotation\_pipeline.py - this program handles read trimming and merging, using FLASH (for merging mate-pairs if analyzing paired-end sequencing) and Trimmomatic for quality control.
* SortMeRNA_example\_script.bash [shell] - example script for using SortMeRNA to remove ribosomal reads from raw data before annotation.  STILL TO COME  

**Annotation:**

* DIAMOND_example\_script.bash - a guide to structuring commmands for DIAMOND to perform annotations or for building dictionaries, including customizable paths to system locations and reference dictionaries.

**Aggregation:**

* DIAMOND_analysis\_counter.py - creates a summary file from DIAMOND results, condensing identical matches and providing counts.  Can be enabled to return either organism or functional results.  STILL TO COME

**Analysis:**    
Note: these programs are located in the "R_scripts" folder.  They all require two sets of DIAMOND summarized results from the aggregation step; an experimental set and a control set.

* diversity_stats.R - this program computes Shannon and Simpson diversity scores.
* make_DESeq\_PCA.R - this program creates a PCA plot of the supplied results (either organism or functional category).
* make_DESeq\_heatmap.R - generates a heatmap contrasting control vs. experimental sets.
* make_combined\_graphs.R - generates two stacked bar graphs: a relative version (y axis stretched to 100%) and an absolute version (y axis is raw counts).
* run_DESeq\_stats.R - computes p-values and most significant differences between summarized results sets.

#### Need assistance?    
Check out the Google Group for SAMSA!  [https://groups.google.com/forum/#!forum/samsa-bioinformatics-group](https://groups.google.com/forum/#!forum/samsa-bioinformatics-group).  

If you find an error with one of these programs, or simply want to ask me questions, you can contact me at [swestreich@gmail.com](mailto:swestreich@gmail.com).  

#### Citations of other tools used:
Westreich, S.T., Korf, I., Mills, D.A., Lemay, D.G.  (2016) SAMSA: A comprehensive metatranscriptome analysis pipeline.  BMC Bioinformatics.

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Magoc, T., and Salzberg, S. FLASH: Fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27:21 (2011), 2957-63.

=======