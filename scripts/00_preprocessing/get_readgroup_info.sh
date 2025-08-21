#!/bin/bash

# change directory
cd /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/fastq/
ls > fastq_files.txt

# create file with list of R1 samples
awk '{print $1}' fastq_files.txt | grep L002_R1_ > R1Samples.txt

# loops through list and print first line
touch sampleReadInfo_snRNA.txt
for sample in `cat R1Samples.txt`; do
   # printf "${sample}\t"
    zcat ${sample} | head -1 >> sampleReadInfo_snRNA.txt	
done;

mv R1Samples.txt  ../scripts/
mv sampleReadInfo_snRNA.txt ../scripts/

cd ../scripts/
paste -d "\t" R1Samples.txt sampleReadInfo_snRNA.txt > sampleReadGroupInfo_snRNA.txt
rm R1Samples.txt
rm sampleReadInfo_snRNA.txt