find -name "*.sra" | awk -F "/" '{print "sam-dump",$0,"|samtools view -bS - >",$2".bam"}' | parallel -j 4 "{}"
