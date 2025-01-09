#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all, kidney, brain and blood sample names
allSamples = list()
type_LBD = list()
LDB_type = list()
sex = list()
read = ["R1", "R2"]
numSamples = 0

with open('sampleReadGroupInfo_snRNA.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace(".", "_")
        line = line.replace("/", "_")
        split = line.split()
        sampleAttributes = split[0].split('_') # 01.RawData/P_65/P_65_CKDL240007500-1A_H5MHCDSXC_S1_L001_R1_001.fastq.gz
        # create a shorter sample name
        stemName = sampleAttributes[2] + "_" + sampleAttributes[3]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "cellranger_dir" : "../cellranger/", 
    "cellbender_dir" : "../cellbender/",
    "cellrangerMAPT_dir" : "../cellrangerMAPT/",
    "fastq_path" : "/research/labs/neurology/fryer/projects/PMI/2024_snRNA/",

    "Comment_Reference" : "This section specifies the location of the mouse reference genome",
    "Mmusculus_dir" : "/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/",
    "Mmusculus_gtf" : "/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf",
    "Mmusculus_fa" : "/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/fasta/genome.fa",
    
    "Comment_MAPT_Reference" : "This section specifies the location of the mouse with human MAPT reference ",
    "MmusculusMAPT_gtf" : "/research/labs/neurology/fryer/projects/references/mouse/human_MAPT_transgenic_mouse_reference/refdata-gex-GRCm39-2024-A_with_human_MAPT/genes/genes.gtf",
    "MmusculusMAPT_fa" : "/research/labs/neurology/fryer/projects/references/mouse/human_MAPT_transgenic_mouse_reference/refdata-gex-GRCm39-2024-A_with_human_MAPT/fasta/genome.fa",
    "MmusculusMAPT_dir" : "/research/labs/neurology/fryer/projects/references/mouse/human_MAPT_transgenic_mouse_reference/refdata-gex-GRCm39-2024-A_with_human_MAPT/cellranger_mouse_ref_with_human_MAPT/",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
'''
outfile.write(header.format(allSamples, read))

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

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace(".", "_")
        line = line.replace("/", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')
        # uniqueNum-number_sequencer_lane_read.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[2] + "_" + sampleAttributes[3]
        stemID = sampleAttributes[2] + '_' + sampleAttributes[3]
        fullName = sampleAttributes[2] + '_' + sampleAttributes[5] + '_' + sampleAttributes[6] + '_' + sampleAttributes[7] 
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00742:821:H5MHCDSXC:1:1101:1597:1000 1:N:0:TGATGATTCA+CGACTCCTAC
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemID  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = fullName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = '10X'

        out = '''
    "{0}":{{
        "fq_path": "/research/labs/neurology/fryer/projects/PMI/2024_snRNA/01.RawData/{7}/",
        "fq1": "{1}",
        "fq2": "{2}",
        "ID": "{3}",
        "SM": "{4}",
        "PU": "{5}",
        "LB": "{6}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, stemName, fullName, PU, LB, stemName))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
