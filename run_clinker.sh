#!/usr/bin/env bash

# conda env
# conda activate ...

# set input names
# WARN: make sure there's no trailing "/"

# INPUTDIR="./input/nife_select1/group4s/overlapping_group4"
# OUTNAME="overlapping_group4_test3"
INPUTDIR="./input/nife_select1/group3s/overlapping_group3"
OUTNAME="overlapping_group3_redo"

# set output names and dirs
name=${INPUTDIR##*/}
outname="./output/${OUTNAME}/${name}"
trimdir="./output/${OUTNAME}/trimmed_gbks"
mkdir -p "./output/${OUTNAME}"
mkdir -p "${trimdir}"

# run clinker with trimming
clinker \
    "$INPUTDIR" \
    -p "${outname}.html" \
    -s "${outname}.json" \
    -mo "${outname}.csv" \
    -o "${outname}_ali.tsv" \
    -dl "," \
    -t "$trimdir" \
    -f

#########################################
# running on trimmed gbk files now...
#########################################

INPUTDIR="${trimdir}"

# run clinker with no trimming
clinker \
    "$INPUTDIR" \
    -p "${outname}.trimmed.html" \
    -s "${outname}.trimmed.json" \
    -mo "${outname}.trimmed.csv" \
    -o "${outname}_ali.trimmed.tsv" \
    -dl "," \
    -f
