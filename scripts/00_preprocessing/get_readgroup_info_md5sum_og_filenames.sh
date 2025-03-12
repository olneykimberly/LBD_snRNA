#!/bin/bash

# change directory
cd /tgen_labs/jfryer/projects/LBD_CWOW/snRNA/tgen_10x/md5sums/
ls > MD5.txt

# create file with list of R1 samples
awk '{print $1}' MD5.txt | grep _R1_ > R1Samples.txt
# remove the .md5 from the file name
sed 's/\.md5$//' R1Samples.txt > R1Samples_snRNA.txt
rm R1Samples.txt

mv R1Samples_snRNA.txt /tgen_labs/jfryer/projects/LBD_CWOW/snRNA/tgen_10x/fastq/R1Samples_snRNA.txt

cd /tgen_labs/jfryer/projects/LBD_CWOW/snRNA/tgen_10x/fastq/
# loops through list and print first line
touch sampleReadInfo_snRNA.txt
for sample in `cat R1Samples_snRNA.txt`; do
   # printf "${sample}\t"
    zcat ${sample} | head -1 >> sampleReadInfo_snRNA.txt	
done;

mv R1Samples_snRNA.txt  /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/R1Samples_snRNA.txt
mv sampleReadInfo_snRNA.txt /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/sampleReadInfo_snRNA.txt

cd /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/
paste -d "\t" R1Samples_snRNA.txt sampleReadInfo_snRNA.txt > sampleReadGroupInfo_snRNA.txt
rm R1Samples_snRNA.txt
rm sampleReadInfo_snRNA.txt

cat sampleReadGroupInfo_snRNA.txt | grep _L002_ > L002_samples_only.txt
rm sampleReadGroupInfo_snRNA.txt 
mv L002_samples_only.txt sampleReadGroupInfo_snRNA.txt
