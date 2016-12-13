# SAMSA, version 2.0

Version 2 of the SAMSA pipeline - faster!  Lighter!  More options!  Less waiting!  Still in beta testing!  

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
(Well, there's nothing here yet.  It's still being built, but stay tuned!)

#### Citations of other tools used:
Westreich, S.T., Korf, I., Mills, D.A., Lemay, D.G.  (2016) SAMSA: A comprehensive metatranscriptome analysis pipeline.  BMC Bioinformatics.

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Magoc, T., and Salzberg, S. FLASH: Fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27:21 (2011), 2957-63.

