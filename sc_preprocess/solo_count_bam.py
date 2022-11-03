# -*- coding: utf-8 -*-
# @Author: Jarning
# @Date:   2021-12-21 19:36:56
# @Last Modified by:   Jarning
# @Last Modified time: 2021-12-27 01:00:11

'''
用途: 用STARsolo对BAM文件直接进行处理，无需转换fastq

1. 目录结构
.
├── bam
│   ├── HumanSpermatogonia_17_1_possorted_genome_bam.bam.1
│   ├── HumanSpermatogonia_17_2_possorted_genome_bam.bam.1
└── samples.txt

2. samples.txt的内容
# BAM文件  样本编号
HumanSpermatogonia_17_1_possorted_genome_bam.bam.1  GSM2928377
HumanSpermatogonia_17_2_possorted_genome_bam.bam.1  GSM2928378
'''

import os
import sys
import re

# 参数设置
CPU = 20
PLATFORM = "10xv2"  # 修改成对应的建库方法，主要用来确定whitelist

inDIR = sys.argv[1]
genome = sys.argv[2]  # GRCh38, mm10
outDIR = inDIR.replace("bam", "solo_outs")

whitelist = {
    "10xv1": "/home/data/data1/resource/singlecell_whitelist/737K-april-2014_rc.txt",
    "10xv2": "/home/data/data1/resource/singlecell_whitelist/737K-august-2016.txt",
    "10xv3": "/home/data/data1/resource/singlecell_whitelist/3M-february-2018.txt.gz"
}

starIndex = {
    "GRCh38": "/home/data/data1/resource/star/GRCh38-2020-A/",
    "mm10": "/home/data/data1/resource/star/mm10-2020-A/"
}

refIndex = starIndex[genome]
WL = whitelist[PLATFORM]

idx = 0
for line in open("samples.txt").readlines():
    idx += 1
    inBAM, gsmID = line.rstrip().split()
    inBAM = os.path.join(inDIR, inBAM)
    outPrefix = os.path.join(outDIR, gsmID)
    solo_cmd = "STAR --genomeDir {refIndex} \
--runThreadN {CPU} \
--outSAMtype BAM Unsorted \
--outFileNamePrefix {outPrefix}/ \
--readFilesIn {inBAM} \
--readFilesCommand samtools view -F 0x100 \
--readFilesType SAM SE \
--soloType CB_UMI_Simple \
--soloCBwhitelist {WL} \
--soloInputSAMattrBarcodeSeq CR UR \
--soloInputSAMattrBarcodeQual CY UY \
--readFilesSAMattrKeep None \
--soloCellFilter EmptyDrops_CR".format(**locals())
    print("# task %s" % idx)
    print(solo_cmd + '\n')
