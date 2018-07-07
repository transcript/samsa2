#!/bin/bash

# Statements common to all scripts
if [[ -z "$SAMSA" ]]; then
  SAMSA="$(readlink -f "${BASH_SOURCE%/*}/..")"
fi
echo "Using SAMSA at $SAMSA" >&2

log()
{
  echo "$@" >&2
}

PROGRAMS="$SAMSA/programs"

PEAR_DIR="$PROGRAMS/pear-0.9.10-linux-x86_64"
PEAR="$PEAR_DIR/bin/pear"
TRIMMOMATIC_DIR="$PROGRAMS/Trimmomatic-0.36"
TRIMMOMATIC="$TRIMMOMATIC_DIR/trimmomatic-0.36.jar"
SORTMERNA_DIR="$PROGRAMS/sortmerna-2.1"
SORTMERNA="$SORTMERNA_DIR/sortmerna"
DIAMOND_DIR="$PROGRAMS"
DIAMOND="$DIAMOND_DIR/diamond"

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
