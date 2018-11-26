#!/bin/bash

# Map the reads against the Atlantic Salmon genome creating:
# - a sam file (GQIxx_salm.sam) with all Salmon reads aligned
# - a fastq file (GQIxx.nosalmreads.fastq) with all non-salmon reads
# created by Joana Meier (Nov 2016) modified by Rishi De-Kayne (July 2017)        

if [[ $1 = "-h" ]]; then
  echo "Usage: Salmon_removal.sh library_raw_reads.fastq GQI_number"
  exit 0
fi

if (( "$#" != 2 )); then
  echo "Please provide 2 arguments: library_raw_reads.fastq GQI_number"
  echo "Usage: Salmon_removal.sh library_raw_reads.fastq"
  exit 0
fi


file=$1
gqi=$2
pu=`awk -F: '{print $1":"$2}' $file | head -1 | sed 's/@//'`

bowtie2 -x /cluster/project/gdc/shared/p298/SalmonGenome/GCF_000233375.1_ICSASG_v2_genomic \
 -q $file --phred33 --end-to-end \
 -p 4 -N 1 --no-unal \
 --un $gqi".noSalmonreads.fastq" \
 -S ""$gqi"_Salmon.sam" \
 --rg-id $gqi"_salm" \
 --rg "ID:"$gqi"_salm" \
 --rg "LB:"$gqi \
 --rg "PL:ILLUMINA" \
 --rg "PU:"$pu \
 --rg "SM:"$gqi"_salm"
