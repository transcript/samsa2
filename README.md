# SAMSA2 - A complete metatranscriptome analysis pipeline

This is a fork of the [original repository](https://github.com/transcript/samsa2). Refer to there for details.

## Modifications and Improvements
The work done on this repo is to adapt the pipeline to work with a Slurm scheduler and some additional parallelisation.
* `bash_scripts/master_scripts.sh` is modified to:
    * use Slurm job arrays
    * utilize tuning parameters in `diamond blastx` (currently tuned to suit 32 threads on 16 CPU cores with hyper threading)
* `python_scripts/DIAMOND_analysis_counter.py` utilizes Python multiprocessing to parallelise the work.

## How To Use
### Preparing input files
The pipeline is designed to work with fastq samples with filenames of the form `*_R{1,2}*.fastq*`. 
`bash_scripts/create_array_dirs.sh` takes a collection of these pairs of fastq files and subdivides them into isolated folders using soft links.
If your files follow a different schema, you may need to modify the script accordingly or create your own.

For example
```{bash}
$ ls -1 input_files
48E_S68_L001_R1_001.fastq
48E_S68_L001_R2_001.fastq
48F_S69_L001_R1_001.fastq
48F_S69_L001_R2_001.fastq
48R1_S73_L001_R1_001.fastq
48R1_S73_L001_R2_001.fastq

$ bash_scripts/create_array_dirs.sh input_files input_files_seperated

$ tree input_files_seperated
./input_files_seperated/
├── 48E_S68_L001
│   ├── 48E_S68_L001_R1_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48E_S68_L001_R1_001.fastq
│   └── 48E_S68_L001_R2_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48E_S68_L001_R2_001.fastq
├── 48F_S69_L001
│   ├── 48F_S69_L001_R1_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48F_S69_L001_R1_001.fastq
│   └── 48F_S69_L001_R2_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48F_S69_L001_R2_001.fastq
├── 48R1_S73_L001
    ├── 48R1_S73_L001_R1_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48R1_S73_L001_R1_001.fastq
    └── 48R1_S73_L001_R2_001.fastq -> /vast/projects/SAMSA/samsa2/input_files/48R1_S73_L001_R2_001.fastq
```

### Modifying input and output dirs
in `bash_scripts/master_script.sh`: 
* modify `INPUT_DIR` to point to the directory with subdivided input files.
* if needed, modify `OUT_DIR` to change the location of where seperate output files are stored.
* if needed, modify `CENTRAL_OUT_DIR` to change the location of where collated output files are stored.

### Submitting the job
Like any other Slurm script, you can submit the master script by
```
sbatch bash_scripts/master_script.sh
```
Note that the most memory intensive step, `diamond blastx`, is currently tuned for performance at the cost of memory. If your nodes do not have the 150G RAM required for the current parameters, you may wish to modify them by reducing the `-b` (block size) and `-c` (number of chunks) values. See the [DIAMOND wiki](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options).
