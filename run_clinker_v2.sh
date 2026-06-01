#!/usr/bin/env bash

# conda env
# conda activate ...
# conda activate ./rp24_scratch2/jamesl2/miniconda/conda/envs/clinker_fork
# pip install -e .

# set input names
# WARN: make sure there's no trailing "/"

# INPUTDIR="./input/nife_select1/group4s/overlapping_group4"
# OUTNAME="overlapping_group4_test3"
# INPUTDIR="./input/nife_select1/group3s/overlapping_group3"
# OUTNAME="overlapping_group3_redo"
# INPUTDIR="./input/01_filtered_sliced_genbanks_pyrospheara/filelist.txt"
# OUTNAME="pyrosphaera_hybrid"

# # set output names and dirs
# name=${INPUTDIR##*/}
outname="./output/${OUTNAME}/${name}"
# trimdir="./output/${OUTNAME}/trimmed_gbks"
# mkdir -p "./output/${OUTNAME}"
# mkdir -p "${trimdir}"

# # run clinker with trimming
# clinker \
#     -fl "$INPUTDIR" \
#     -p "${outname}.html" \
#     -s "${outname}.json" \
#     -mo "${outname}.csv" \
#     -o "${outname}_ali.tsv" \
#     --svg "${outname}.svg" \
#     --anchor_target_genes \
#     -dl "," \
#     -t "$trimdir" \
#     -f

#########################################
# running on trimmed gbk files now...
#########################################

OUTNAME="pyrosphaera_hybrid"
outname="./output/${OUTNAME}"
trimdir="./output/pyro_trimmed_gbks"
TRIMINPUT="${trimdir}/filelist.txt"

# generate file list
fd . "$trimdir" -e gbk -a >$TRIMINPUT

# run clinker with no trimming
clinker \
    -fl "$TRIMINPUT" \
    -p "${outname}.trimmed.html" \
    -s "${outname}.trimmed.json" \
    -mo "${outname}.trimmed.csv" \
    -o "${outname}_ali.trimmed.tsv" \
    --svg "${outname}.trimmed.svg" \
    --anchor_target_genes \
    -dl "," \
    -f
