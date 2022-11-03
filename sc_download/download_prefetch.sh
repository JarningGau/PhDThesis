cat todo.txt | parallel -j 4 "prefetch {} -X 200G -O ./"
