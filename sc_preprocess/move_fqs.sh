find -name "*.gz" | awk -F '/' '{print "mv",$0,$2"/"$4}' | bash
find -name "*Missing*" | xargs rmdir

cat meta.tsv | awk '{print $2}'|sort|uniq|xargs mkdir
cat meta.tsv | awk '{print "mv",$1"_1.fastq.gz",$2}'|bash
cat meta.tsv | awk '{print "mv",$1"_2.fastq.gz",$2}'|bash

