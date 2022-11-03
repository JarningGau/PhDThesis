# -*- coding: utf-8 -*-
# @Author: Jarning
# @Date:   2021-12-21 19:36:56
# @Last Modified by:   Jarning
# @Last Modified time: 2022-07-18 15:45:51

'''
用途: 用STARsolo对FASTQ文件进行预处理

1. 目录结构
.
└── fastq
    └── GSM3526583
        ├── SRR8363216_1.fastq.gz
        ├── SRR8363216_2.fastq.gz
        ├── SRR8363217_1.fastq.gz
        └── SRR8363217_2.fastq.gz
2. demo
# INPUT: path contains sample dirs
# OUTPUT: solo count bash commands
python solo_count.fastq.py /path/to/fastq > solo_count.sh


'''

import os
import sys
import re

# 参数设置
CPU = 20
PLATFORM = "10xv2"  # 10xv2, 10xv3, dropseq, STRTseq_Tang

inDIR = sys.argv[1]
genome = "GRCh38"  # GRCh38, mm10, mmul10
outDIR = inDIR.replace("fastq", "solo_outs")
R1 = re.compile(r'.+_R{0,1}1_{0,1}.*?.fastq.gz$')
R2 = re.compile(r'.+_R{0,1}2_{0,1}.*?.fastq.gz$')
R3 = re.compile(r'.+_R{0,1}3_{0,1}.*?.fastq.gz$')

# outSAMSettings = "--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM" # output bam
outSAMSettings = "--outSAMtype None"  # no bam output

whitelist = {
    "10xv1": "/home/data/data1/resource/singlecell_whitelist/737K-april-2014_rc.txt",
    "10xv2": "/home/data/data1/resource/singlecell_whitelist/737K-august-2016.txt",
    "10xv3": "/home/data/data1/resource/singlecell_whitelist/3M-february-2018.txt",
    "dropseq": "None",
    "STRTseq_Tang": "/home/data/data1/resource/singlecell_whitelist/Tang_STRTseq_96_barcodes.txt"
}

# for STATseq, see https://github.com/JarningGau/NGS_pipeline/tree/master/scrna_pipeline CB+UMI(8+8) on reads2
# for dropseq, see https://github.com/alexdobin/STAR/issues/584 CB+UMI(12+8) on reads1
# for 10x, If input 10x FQ is not standard format (PE100 or PE150), add --soloBarcodeReadLength 0.
barcode_setting = {
    "10xv2": "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10",
    "10xv3": "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12",
    "dropseq": "--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --soloBarcodeReadLength 0",
    "STRTseq_Tang": "--soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 8 --soloBarcodeReadLength 0"
}

starIndex = {
    "GRCh38": "/home/data/data1/resource/star/GRCh38-2020-A/",
    "mm10": "/home/data/data1/resource/star/mm10-2020-A/",
    "mmul10": "/home/data/data1/resource/star/Mmul_10"
}

refIndex = starIndex[genome]
BS = barcode_setting[PLATFORM]
WL = whitelist[PLATFORM]


def get_paired_fq(path):
    inFASTQ_cDNA, inFASTQ_barcode = '', ''
    fq1 = sorted([os.path.join(path, f) for f in os.listdir(path) if R1.match(f)])
    fq2 = sorted([os.path.join(path, f) for f in os.listdir(path) if R2.match(f)])
    fq3 = sorted([os.path.join(path, f) for f in os.listdir(path) if R3.match(f)])
    if len(fq1) != len(fq2):
        raise SystemExit("The PE fastq files were not matched!")
    if len(fq1) == 0 or len(fq2) == 0:
        raise SystemExit("No fastq file was found. Please check the file name.")
    if len(fq3) == 0:
        if PLATFORM in ["10xv2", "10xv3", "dropseq"]:
            inFASTQ_cDNA = ",".join(fq2)
            inFASTQ_barcode = ",".join(fq1)
        elif PLATFORM in ["STRTseq_Tang"]:
            inFASTQ_cDNA = ','.join(fq1)
            inFASTQ_barcode = ",".join(fq2)
        else:
            raise SystemExit('%s was not support! Please choose correct library!' % PLATFORM)
    elif len(fq3) == len(fq2) == len(fq1):
        if PLATFORM in ["10xv2", "10xv3", "dropseq"]:
            inFASTQ_cDNA = ','.join(fq3)
            inFASTQ_barcode = ",".join(fq2)
    else:
        raise SystemExit("The fastq files was not properly matched!")
    return inFASTQ_cDNA, inFASTQ_barcode


idx = 0
for sample in sorted(os.listdir(inDIR)):
    idx += 1
    inFASTQ_cDNA, inFASTQ_barcode = get_paired_fq(os.path.join(inDIR, sample))
    outPrefix = os.path.join(outDIR, sample)
    solo_cmd = "STAR --genomeDir {refIndex} \
--runThreadN {CPU} \
{outSAMSettings} \
--outFileNamePrefix {outPrefix}/ \
--readFilesIn {inFASTQ_cDNA} {inFASTQ_barcode} \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
{BS} \
--soloCBwhitelist {WL} \
--soloCellFilter EmptyDrops_CR".format(**locals())
    print("#### task %s" % idx)
    print(solo_cmd + '\n')
