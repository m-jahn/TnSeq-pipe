#!/usr/bin/env bash

# input parameters: reference genome, primer model
ref=${ref:-"ASM928v2"}
model=${model:-model_pKMW7}
data=${data:-"data/fastq/"}
pattern=${pattern:-".*"}
output=${output:-"data/"}
stepSize=${stepSize:-11}
tileSize=${tileSize:-11}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
    fi
    shift
done

# assign paths. The subdirectories to data in and output are
# currently hardcoded
MAPTN="feba/bin/MapTnSeq.pl"
TNPOOL="feba/bin/DesignRandomPool.pl"
REF="ref/"
MODEL="feba/primers/$model"
FASTQ=$data
OUT=$output"mapped/"
POOL=$output"pool/"

# if output folders are not present, create them
for dir in $OUT $POOL
do
    if [ -d ${dir} ]; then
        echo "Output directory: ${dir} exists"
    else
        echo "Output directory: ${dir} does not exist, creating.."
        mkdir ${dir:0:-1};
    fi
done

# retrieve genome files from NCBI
Rscript source/prepare_ref_genome.R ${ref}
fasta_path=`ls ${REF} | grep ${ref}.*.fna`
gff_path=`ls ${REF} | grep ${ref}.*.gff`
[ -f ${REF}/${fasta_path} ] || echo "FASTA file was not found in the specified path: ${REF}/${fasta_path}"
[ -f ${REF}/${gff_path} ] || echo "GFF file was not found in the specified path: ${REF}/${gff_path}"

# make filename pattern and print file number
filepattern=${pattern}.fastq.gz$
echo "Input file pattern matches: ${filepattern}"
filenames=(`ls ${FASTQ} | grep ${filepattern}`)
echo "Input files matching pattern: ${#filenames[@]}"
if [ ${#filenames[@]} == 0 ]; then
    echo "Found no files to process. Exiting"
    exit
fi


# step 1: trim reads and map barcodes to reference genome
# read file names
ls ${FASTQ} | grep ${filepattern} | while read fastq;
do
    # extract ID of fastq.gz file
    ID=`echo ${fastq} | cut -f 1 -d \.`
    
    perl ${MAPTN} \
    -stepSize ${stepSize} -tileSize ${tileSize} \
    -genome ${REF}/${fasta_path} \
    -model ${MODEL} \
    -first ${FASTQ}/${fastq} \
    -unmapped ${OUT}/${ID}_unmapped.txt \
    -trunc ${OUT}/${ID}_truncated.txt \
    > ${OUT}${ID}.tsv
done  | parallel --no-notice --bar
# optional piping to parallel jobs.
# use --jobs N to set max number of CPU cores


# step 2: create summary table, include all barcodes (min  = 1)
# can filter by N barcdes later
perl ${TNPOOL} \
-minN 1 \
-pool ${POOL}pool.tsv \
-genes ${REF}/${gff_path} \
${OUT}*.tsv
