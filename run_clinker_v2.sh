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
INPUTLIST="/home/james/Dropbox/Monash/HydDB/famous_hydrogenases/Group2a_NiFes/01_group2a_filtered_sliced_genbanks/spacedust/output/plot/huc_clustersearch_clustsize3_cov70-filelist.txt"
OUTNAME="group2a_spacedust_to_clinker"
# INPUTLIST="/home/james/Dropbox/Monash/Pyrosphaera/gene_neighbourhood_extraction/Pyrosphaera_fefenife_hybrid_neighbours/02_filtered_and_sliced_gbks_5000/spacedust/output/plot/pyro_clustersearch.clustsize3.cov70.minseqid20-clusterIDs.txt"
# OUTNAME="pyro_spacedust_to_clinker_clustsize3"

# # set output names and dirs
name=${INPUTLIST##*/}
outname="./output/${OUTNAME}/${name}"
trimdir="./output/${OUTNAME}/trimmed_gbks"
mkdir -p "./output/${OUTNAME}"
mkdir -p "${trimdir}"

# run clinker with trimming
clinker \
    -fl "$INPUTLIST" \
    -p "${outname}.html" \
    -s "${outname}.json" \
    -mo "${outname}.csv" \
    -o "${outname}_ali.tsv" \
    --svg "${outname}.svg" \
    --anchor_target_genes \
    --identity 0.5 \
    -dl "," \
    -t "$trimdir" \
    -na \
    -f

# use '-na' for no alignment... should add colour_map though

#########################################
# running on trimmed gbk files now...
#########################################

# OUTNAME="pyrosphaera_hybrid"
# outname="./output/${OUTNAME}"
# trimdir="./output/pyro_trimmed_gbks"
# TRIMINPUT="${trimdir}/filelist.txt"
#
# # generate file list
# fd . "$trimdir" -e gbk -a >$TRIMINPUT
#
# # run clinker with no trimming
# clinker \
#     -fl "$TRIMINPUT" \
#     -p "${outname}.trimmed.html" \
#     -s "${outname}.trimmed.json" \
#     -mo "${outname}.trimmed.csv" \
#     -o "${outname}_ali.trimmed.tsv" \
#     --svg "${outname}.trimmed.svg" \
#     --anchor_target_genes \
#     -dl "," \
#     -f
