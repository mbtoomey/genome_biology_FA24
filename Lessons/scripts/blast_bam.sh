samtools view -b -S -T CYP2J19.fasta HEK_CYP2J19_blast.sam -o HEK_CYP2J19_blast.bam

samtools sort HEK_CYP2J19_blast.bam -o HEK_CYP2J19_blast_sort.bam

samtools index HEK_CYP2J19_blast_sort.bam