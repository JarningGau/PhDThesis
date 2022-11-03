find -name "*.sra" | parallel -j 8 "fastq-dump --gzip --split-files {}"
