#!/bin/bash

# Create variables for commonly-used utilities; allow environment override
MKDIR=${MKDIR:-'mkdir -p'}
RMDIR=${RMDIR:-rmdir}
CD=${CD:-cd}
RM=${RM:-rm}
MV=${MV:-mv}
TOUCH=${TOUCH:-touch}
PYTHON=${PYTHON:-python}
JAVA=${JAVA:-java}
GUNZIP=${GUNZIP:-gunzip}

# Statements common to all scripts
if [[ -z "$SAMSA" ]]; then
  SAMSA="$(readlink -f "${BASH_SOURCE%/*}/..")"
fi
echo "Using SAMSA at $SAMSA" >&2

# Functions used by various scripts
log() {
  echo "$@" >&2
}

checked() {
  $@
  status=$?
  if [[ $status -ne 0 ]]; then
    echo "'$@' exited with non-zero status $status" >&2
    exit $status
  fi
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
export R_LIBS="$SAMSA/R_scripts/packages"

# Unless indicated otherwise, ensure the utilities are usable
if [[ -z "$IGNORE_PATHS" ]]; then
  if [[ ! -f "$PEAR" ]]; then
    echo "ERROR: PEAR not found (did you extract it?) at $PEAR" >&2
    exit 1
  fi

  if [[ ! -f "$TRIMMOMATIC" ]]; then
    echo "ERROR: Trimmomatic not found (did you extract it?) at $TRIMMOMATIC" >&2
    exit 1
  fi

  if [[ ! -f "$SORTMERNA" ]]; then
    echo "ERROR: SortMeRNA not found (did you extract and build it?) at $SORTMERNA" >&2
    exit 1
  fi

  if [[ ! -f "$DIAMOND" ]]; then
    echo "ERROR: Diamond not found (did you extract it?) at $DIAMOND" >&2
    exit 1
  fi
fi
