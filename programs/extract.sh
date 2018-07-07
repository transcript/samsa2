#!/bin/bash

# Script to extract the programs in this directory

tar xzf "$(dirname $0)/diamond-linux64.tar.gz"
tar xzf "$(dirname $0)/pear-0.9.10-linux-x86_64.tar.gz"
tar xzf "$(dirname $0)/sortmerna-2.1.tar.gz"
unzip "$(dirname $0)/Trimmomatic-0.36.zip"
