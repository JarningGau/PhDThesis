Copy these files to a dir to store your download files.

meta.tsv             # from SRA database
# cat meta.tsv
# SRR6342093	GSM2874055
# SRR6342094	GSM2874056
SRR_Acc_List.txt     # from SRA database
download_prefetch.sh
sra2fq.sh

1. Dowload data form SRA database.
cat SRR_Acc_List.txt > todo.txt # prepare todo list
bash download_prefetch.sh       # download via sra toolkit in parallel
# this step will generate .sra files under current path.
#├── SRR6342193
#│   └── SRR6342193.sra
#├── SRR6342194
#│   └── SRR6342194.sra
bash sra2fq.sh                  # .sra -> . fastq.gz
# this step will generate *.fastq.gz files under current path
#├── SRR6342193_1.fastq.gz
#├── SRR6342193_2.fastq.gz
#├── SRR6342194_1.fastq.gz
#├── SRR6342194_2.fastq.gz

2. Rename the *.fastq.gz files
# create sample dir (GSM_ID) and move *.fastq.gz files to it
cat meta.tsv | awk '{print $2}'|sort|uniq|xargs mkdir     # create sample dir (GSM_ID)
cat meta.tsv | awk '{print "mv",$1"_1.fastq.gz",$2}'|bash # move _1.fastq.gz
cat meta.tsv | awk '{print "mv",$1"_2.fastq.gz",$2}'|bash # move _2.fastq.gz
