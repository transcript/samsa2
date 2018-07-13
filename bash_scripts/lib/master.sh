#!/bin/bash

#######################################################################
#
# lib/master.sh
# Modular master_script component definitions.
#
# This script encapsulates each pipeline step in a sequence of bash
# functions for use in both the master_script and its test script.
#
# Expected usage for direct invocation:
# SAMSA_INTERACTIVE=1 source bash_scripts/lib/master.sh
#
#######################################################################

if [[ -z "$SAMSA" ]]; then
  echo "WARNING: SAMSA environment not initialized" >&2
  source "${BASH_SOURCE%/*}/../bash_scripts/lib/common.sh"
fi

INPUT_DIR="${INPUT_DIR:-$SAMSA/input_files}"

DIA_REFSEQ="${DIA_REFSEQ:-$SAMSA/full_databases/RefSeq_bac}"
DIA_SUBSYS="${DIA_SUBSYS:-$SAMSA/full_databases/subsys_db}"
DIA_REFSEQ_PATH="${DIA_REFSEQ}.fa"
DIA_SUBSYS_PATH="${DIA_SUBSYS}.fa"

# do_pear <input-path> <output-path>
do_pear() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  for r1f in $ipath/*_R1*; do
    r2f="$(echo $r1f | sed 's/_R1/_R2/')"
    obase="$(repath "$r1f" "$opath" "s|_R1.*||").merged"
    checked $PEAR -f $r1f -r $r2f -o $obase
    # output: ${obase}.merged.fastq
  done
}

# do_trimmomatic <input-path> <output-path>
do_trimmomatic() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  for mf in $ipath/*.merged.*; do
    ofile="$(repath "$mf" "$opath" 's|\.merged|.cleaned.fastq|')"
    do_java -jar $TRIMMOMATIC SE -phred33 $mf $ofile SLIDINGWINDOW:4:15 MINLEN:99
  done
}

# do_read_counter <input-path> <output-file>
do_read_counter() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  if [[ -f "$opath/read_counts.txt" ]]; then
    rm "$opath/read_counts.txt";
  fi
  touch "$opath/read_counts.txt"
  for f in $ipath/*.cleaned.fastq; do
    do_python $PY_DIR/raw_read_counter.py -I "$f" -O "$opath/read_counts.txt"
  done
}

# do_sortmerna <input-path> <output-path>
do_sortmerna() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  local fasta_db="$SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta"
  local index_db="$SORTMERNA_DIR/index/silva-bac-16s-db"
  for f in $ipath/*.cleaned.fastq; do
    fbase="$(echo $f | sed 's|\.fastq$||')"
    obase="$(repath "$fbase" "$opath" 's|\.cleaned|.ribodepleted|')"
    checked $SORTMERNA \
      --ref "$fasta_db,$index_db" \
      --reads $f \
      --aligned ${fbase}.ribosomes \
      --other $obase \
      --fastx \
      --num_alignments 0 \
      --log -v
  done
}

# do_diamond_refseq <input-path> <output-path>
do_diamond_refseq() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  do_mkdir "$opath/daa_binary_files"
  for f in $ipath/*.ribodepleted.fastq; do
    dbase="$(repath "$f" "$opath/daa_binary_files" 's|ribodepleted\.fastq|RefSeq|')"
    ofile="$(repath "$f" "$opath" 's|ribodepleted\.fastq|RefSeq_annotated|')"
    checked $DIAMOND blastx --db "$DIA_REFSEQ" -q "$f" -a "$dbase" -t . -k 1
    checked $DIAMOND view --daa "${dbase}.daa" -o "$ofile" -f tab
  done
}

# do_diamond_subsys <input-path> <output-path>
do_diamond_subsys() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  do_mkdir "$opath/daa_binary_files"
  for f in $ipath/*.ribodepleted.fastq; do
    dbase="$(repath "$f" "$opath/daa_binary_files" 's|ribodepleted\.fastq|Subsys|')"
    ofile="$(repath "$f" "$opath" 's|ribodepleted\.fastq|subsys_annotated|')"
    checked $DIAMOND blastx --db "$DIA_SUBSYS" -q "$f" -a "$dbase" -t . -k 1
    checked $DIAMOND view --daa "${dbase}.daa" -o "$ofile" -f tab
  done
}

# do_refseq_analysis <input-path> <output-path>
do_refseq_analysis() {
  local ipath="$1"
  local opath="$2"
  for f in $ipath/*RefSeq_annotated*; do
    # FIXME: Add output directory argument to DIAMOND_analysis_counter.py
    do_python "$PY_DIR/DIAMOND_analysis_counter.py" -I "$f" -D "$DIA_REFSEQ_PATH" -O
    do_python "$PY_DIR/DIAMOND_analysis_counter.py" -I "$f" -D "$DIA_REFSEQ_PATH" -F
  done
  do_mkdir "$opath/RefSeq_results/org_results"
  do_mkdir "$opath/RefSeq_results/func_results"
  do_mv $ipath/*organism.tsv $opath/RefSeq_results/org_results
  do_mv $ipath/*function.tsv $opath/RefSeq_results/func_results
}

# do_subsys_analysis <input-path> <output-path>
do_subsys_analysis() {
  local ipath="$1"
  local opath="$2"
  do_mkdir "$opath"
  do_mkdir "$opath/Subsystems_results/receipts"
  for f in $ipath/*subsys_annotated*; do
    ofile="$(repath "$f" "$opath" 's|ribodepleted\.fastq|RefSeq_annotated|')"
    do_python $PY_DIR/DIAMOND_subsystems_analysis_counter.py \
      -I "$f" \
      -D "$DIA_SUBSYS_PATH" \
      -O "${ofile}.hierarchy" \
      -P "${ofile}.receipt"
    do_python $PY_DIR/subsys_reducer.py \
      -I "${ofile}.hierarchy"
  done
  do_mv $opath/*.reduced "$opath/Subsystems_results"
  do_mv $opath/*.receipt "$opath/Subsystems_results/receipts"
  # NOTE: Artifacts $opath/*.hierarchy are kept
}

# do_deseq <input-path> <output-path> <raw_counts-path>
do_deseq() {
  local ipath="$1"
  local opath="$2"
  local cpath="$3"
  do_mkdir "$opath"
  do_rscript $R_DIR/run_DESeq_stats.R \
    -I "$ipath/RefSeq_results/org_results" \
    -O "$opath/RefSeq_org_DESeq_results.tab" \
    -R "$cpath"
  do_rscript $R_DIR/run_DESeq_stats.R \
    -I "$ipath/RefSeq_results/func_results" \
    -O "$opath/RefSeq_func_DESeq_results.tab" \
    -R "$cpath"
  do_rscript $R_DIR/Subsystems_DESeq_stats.R \
    -I "$ipath/Subsystems_results" \
    -O "$opath/Subsystems_level-1_DESeq_results.tab" \
    -R "$cpath" \
    -L 1
}

