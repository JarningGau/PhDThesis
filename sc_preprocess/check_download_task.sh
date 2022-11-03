# done sra
ls | grep .sra | awk -F '.' '{print $1}' > done.txt
# check todo task
cat SRR_Acc_List.txt done.txt | sort | uniq -c | grep "1 " | awk '{print $2}'
# group and move fastq
cat meta.tsv |awk '{print "mv",$1"_2.fastq.gz",$2}'
