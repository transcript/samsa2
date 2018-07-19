#!/bin/bash

# Common library script used by shell scripts in SAMSA2
#
# Expected usage for scripts in nested directories:
# source "${BASH_SOURCE%/*}/../bash_scripts/lib/common.sh"
#
# Expected usage for direct invocation:
# SAMSA_INTERACTIVE=1 source bash_scripts/lib/common.sh
#
# Environment configuration:
#   VARIABLE      DESCRIPTION
#   ------------- ---------------------------------------------------------
#   SAMSA         Top-level SAMSA2 directory override
#   DEBUG         When non-empty, output debugging and progress information
#   LOGFILE       Logging file path override
#   IGNORE_DEPS   If set, ignore the dependency check for included packages
#   DRY_RUN       If set, just echo what would be done

# Create variables for commonly-used utilities; allow environment override
MKDIR=${MKDIR:-'mkdir -p'}
RMDIR=${RMDIR:-rmdir}
RM=${RM:-rm}
MV=${MV:-mv}
TOUCH=${TOUCH:-touch}
PYTHON=${PYTHON:-python}
JAVA=${JAVA:-java}
R=${R:-R}
RSCRIPT=${RSCRIPT:-Rscript}
TAR=${TAR:-tar}
GUNZIP=${GUNZIP:-gunzip}

# warn <msg...>
warn() {
  logstr "WARNING: $@"
}

# error <msg...>
error() {
  logstr "ERROR: $@"
}

# fatal <msg...>
fatal() {
  error "$@"
  if [[ -z "$SAMSA_INTERACTIVE" ]]; then
    exit 1
  fi
}

# log <cmd...>
log() {
  "$@" 2>&1 | tee -a "$LOGFILE"
  return ${PIPESTATUS[0]}
}

# logstr <msg...>
logstr() {
  echo "$@" | tee -a "$LOGFILE"
}

# debug <msg...>
debug() {
  if [[ -n "$DEBUG" ]]; then
    logstr "+ $(date +%H:%m:%S) $@"
  fi
}

# checked <cmd...>
checked() {
  if [[ -n "$DRY_RUN" ]]; then
    debug "DRY RUN: $@"
    logstr "$@"
  else
    debug "Running $@"
    log "$@"
    status=$?
    if [[ $status -ne 0 ]]; then
      echo "'$@' exited with non-zero status $status" >&2
      exit $status
    fi
    debug "Finished running $@"
  fi
}

# chk_value <value> <emsg>
chk_value() {
  if [[ -z "$1" ]]; then
    error "$2: no value"
    return 1
  fi
}

# chk_file <path> <emsg>
chk_file() {
  if [[ ! -f "$1" ]]; then
    error "$2: $1"
    return 1
  fi
}

# chk_dir <path> <emsg>
chk_dir() {
  if [[ ! -d "$1" ]]; then
    error "$2: $1"
    return 1
  fi
}

do_mkdir() { checked $MKDIR $@; }
do_rmdir() { checked $RMDIR $@; }
do_rm() { checked $RM $@; }
do_mv() { checked $MV $@; }
do_touch() { checked $TOUCH $@; }
do_python() { checked $PYTHON $@; }
do_java() { checked $JAVA $@; }
do_r() { checked $R $@; }
do_rscript() { checked $RSCRIPT $@; }
do_tar() { checked $TAR $@; }
do_gunzip() { checked $GUNZIP $@; }

# Path to top-level directory
if [[ -z "$SAMSA" ]]; then
  SAMSA="$(readlink -f "${BASH_SOURCE%/*}/../..")"
fi
echo "Using SAMSA at $SAMSA" >&2

# Logging
if [[ ! -d "$SAMSA/logs" ]]; then
  $MKDIR "$SAMSA/logs"
fi
if [[ -z "$LOGFILE" ]]; then
  LOGFMT="$SAMSA/logs/out-$(date +%Y%m%d%H%m%S)-%d.log"
  LOGIDX=0
  LOGFILE="$(printf $LOGFMT $LOGIDX)"
  while [[ -f "$LOGFILE" ]]; do
    LOGIDX=$(($LOGIDX + 1))
    LOGFILE="$(printf $LOGFMT $LOGIDX)"
  done
fi

export LOGFILE

# Functions used by various scripts

# repath <path> <new-dir> [regex]
repath() {
  local ifile="$1"
  local opath="$2"
  local re="$3"
  local fname="$(basename "$ifile")"
  if [[ -n "$re" ]]; then
    fname="$(echo $fname | sed "$re")"
  fi
  echo $fname
}

# Paths to included components
PROGRAMS="$SAMSA/programs"
PY_DIR="$SAMSA/python_scripts"
R_DIR="$SAMSA/R_scripts"

# Paths to included utilities
PEAR_DIR="$PROGRAMS/pear-0.9.10-linux-x86_64"
PEAR="$PEAR_DIR/bin/pear"
TRIMMOMATIC_DIR="$PROGRAMS/Trimmomatic-0.36"
TRIMMOMATIC="$TRIMMOMATIC_DIR/trimmomatic-0.36.jar"
SORTMERNA_DIR="$PROGRAMS/sortmerna-2.1"
SORTMERNA="$SORTMERNA_DIR/sortmerna"
DIAMOND_DIR="$PROGRAMS"
DIAMOND="$DIAMOND_DIR/diamond"

# Ensure locally-downloaded R packages are usable
if [[ -n "$R_LIBS" ]]; then
  R_LIBS="$R_LIBS:$SAMSA/R_scripts/packages"
else
  R_LIBS="$SAMSA/R_scripts/packages"
fi
export R_LIBS

# Unless indicated otherwise, ensure the utilities are usable
if [[ -z "$IGNORE_DEPS" ]]; then
  if [[ ! -f "$PEAR" ]]; then
    fatal "PEAR not found (did you extract it?) at $PEAR"
  fi

  if [[ ! -f "$TRIMMOMATIC" ]]; then
    fatal "Trimmomatic not found (did you extract it?) at $TRIMMOMATIC"
  fi

  if [[ ! -f "$SORTMERNA" ]]; then
    fatal "SortMeRNA not found (did you extract and build it?) at $SORTMERNA"
  fi

  if [[ ! -f "$DIAMOND" ]]; then
    fatal "Diamond not found (did you extract it?) at $DIAMOND"
  fi
fi

