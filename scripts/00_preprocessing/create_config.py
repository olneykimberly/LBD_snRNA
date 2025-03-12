#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all sample names
allSamples = list()
numSamples = 0

with open('sampleReadGroupInfo_snRNA.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        #line = line.replace(".", "_")
        #line = line.replace("/", "_")
        split = line.split()
        sampleAttributes = split[0].split('_') # MAYO_0407_1_BR_Nuclei_C1_X3SC3_A16845_22KJLTLT4_AGCAAGAAGC_L007_R1_001.fastq.gz
        # [0] study
        # [1] patient
        # [2] visit (1 for frozen banked samples)
        # [3] source (2 letter code; BR = brain)
        # [4] fraction (Nuclie or whole)
        # [5] subgroup - C1 is control as they don't know the disease status
        # [6] assay 
        # [7] library (1 letter and 5 numbers)
        # [8] sequence order - was [8:9] barcode 
        # [9] lane
        # [10] strand or index
      
        # create a shorter sample name
        stemName = sampleAttributes[3] + "_" + sampleAttributes[4] + "_" + sampleAttributes[1]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "cellranger_dir" : "../cellranger/", 
    "cellbender_dir" : "../cellbender/",
    "fastq_path" : "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/fastq/",

    "Comment_Reference" : "This section specifies the location of the reference genome",
    "human_dir" : "/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A/",
    "human_gtf" : "/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A/genes/genes.gtf",
    "human_fa" : "/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa",
    
    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
'''
outfile.write(header.format(allSamples))

# config formatting
counter = 0
with open('sampleReadGroupInfo_snRNA.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        base = base.replace("_R1_", "")
        sampleInfo = split[1]

        split = line.split()
        sampleAttributes = split[0].split('_')
        # uniqueNum-number_sequencer_lane_read.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[3] + "_" + sampleAttributes[4] + "_" + sampleAttributes[1]
        seq_order = sampleAttributes[0] + "_" + sampleAttributes[1] + "_" + sampleAttributes[2] + "_" + sampleAttributes[3] + "_" + sampleAttributes[4] + "_" + sampleAttributes[5] + "_" + sampleAttributes[6] + "_" + sampleAttributes[7] + "_" + sampleAttributes[8]
        stemID = sampleAttributes[3] + "_" + sampleAttributes[4] + "_" + sampleAttributes[1]
        fullName = sampleAttributes[0] + "_" + sampleAttributes[1] + "_" + sampleAttributes[2] + "_" + sampleAttributes[3] + "_" + sampleAttributes[4] + "_" + sampleAttributes[5] + "_" + sampleAttributes[6] + "_" + sampleAttributes[7] 
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @LH00295:158:22KJLTLT4:2:1101:3743:1042 1:N:0:GTGGATCAAA+CAGGGTTGGC
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        SEQ = seq_order
        ID = stemID  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = fullName  # Sample
        PU = flowcellID  # Platform Unit
        LB = '10X'

        out = '''
    "{0}":{{
        "fq_path": "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/fastq/",
        "SEQ": "{1}",
        "ID": "{2}",
        "SM": "{3}",
        "PU": "{4}",
        "LB": "{5}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, seq_order, stemName, fullName, PU, LB, stemName))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
