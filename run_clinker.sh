#!/usr/bin/env bash

# conda env
# conda activate ...

# set input names
# WARN: make sure there's no trailing "/"

INPUTDIR="./input/group2a_sliced_genbanks"
OUTNAME="group2a_sliced_round1"

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

#####

INPUTDIR=$OUTNAME
OUTNAME="group2a_sliced_round2"

# set output names and dirs
name=${INPUTDIR##*/}
outname="./output/${OUTNAME}/${name}"
trimdir="./output/${OUTNAME}/trimmed_gbks"
mkdir -p "./output/${OUTNAME}"
mkdir -p "${trimdir}"

# run clinker with no trimming
clinker \
    "$INPUTDIR" \
    -p "${outname}.html" \
    -s "${outname}.json" \
    -mo "${outname}.csv" \
    -o "${outname}_ali.tsv" \
    -dl "," \
    -f
