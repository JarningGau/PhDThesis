# -*- coding: utf-8 -*-
# @Author: Jarning
# @Date:   2021-12-20 17:34:33
# @Last Modified by:   Jarning
# @Last Modified time: 2021-12-20 17:55:19

import os
import sys
import re

inDIR = sys.argv[1]
outDIR = inDIR.replace("fastq", "kb_outs")
R1 = re.compile(r'.+_R1_.+.fastq.gz')
R2 = re.compile(r'.+_R2_.+.fastq.gz')

CPU = 20
MEM = "48G"
PLATFORM = "10xv2"
refIndex = "/home/data/data1/resource/kb/GRCh38-2020-A/index.idx"
refT2G = "/home/data/data1/resource/kb/GRCh38-2020-A/t2g.txt"


def get_paired_fq(path):
    fq1 = sorted([os.path.join(path, f) for f in os.listdir(path) if R1.match(f)])
    fq2 = sorted([os.path.join(path, f) for f in os.listdir(path) if R2.match(f)])
    return ' '.join([a + " " + b for a, b in zip(fq1, fq2)])


FASTQs = get_paired_fq(inDIR)
kb_cmd = "kb count -i {refIndex} -g {refT2G} -t {CPU} -m {MEM} -x {PLATFORM} --tmp ./tmp --filter bustools --cellranger -o {outDIR} {FASTQs}".format(
    **locals())
print(kb_cmd)
