# RNA-Seq analysis

The sequencing of RNA transcripts has become the standard method for untargeted quantification of gene expression. Unlike genome sequencing, where we aim to have uniform coverage or sequencing depth across the entire geneome, in RNA-seq we expect both the coverage and depth of sequencing to vary in proportion to pattern and levels of expression of each transcript in each sample. Thus, we can quantify gene expression levels by counting the number of sequencing reads that map to a reference transcription or genome. 

RNA-Seq sample preparation differs from DNA sequencing preparation in several ways: 
- RNA is more unstable than DNA, requiring additional care during sampling and processing
- There are many forms of RNA in the cell and much of it has functions other than protein expression (e.r. ribosomal RNAs). Therefore, mosgt experiments use some for of selection to enrich for specific types of RNAs
- RNA must be reverse transcribed to cDNA before short-read sequencing libraries can be assembled. 

![](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863231/bin/nihms768779f1.jpg)
Diagram from: 10.1101/pdb.top084970

Once libraries are made, they can be sequenced in much the same ways as DNA libraries. The resulting sequence reads can then be aligned to a reference transcriptome or annotated genome, read counts generated, and then exported for statistical analyses of expression patterns. 

![](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863231/bin/nihms768779f2.jpg)

## The Experiment

For this example we will analyze RNA-Seq data from an experiment in my lab. Recently, [we have discovered](https://www.sciencedirect.com/science/article/pii/S0960982222012908) that the enzymes CYP2J19 and BDH1L catalyzes the conversion of common yellow carotenoid pigments in the diets of birds into red ketocarotenoids that color their feathers and skin.

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image1.png)

We have also discovered that the activity of these enzymes in enhanced by the expression of the protein TTC39B.

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image2.png)

However, we do not know how TTC39B is enhancing carotenoid metabolism. Therefore, we wanted to know if and how **TTC39B** expression alters gene expression within cells. To do this we cultured [HEK293 cells] and transfected these cells with expression constructs for **CYP2J19** and **BDH1L** or **CYP2J19**, **BDH1L**, and **TTC39B**

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image3.png)

After 48 hours of expression, we harvested these cells, extracted RNA, and prepare short-read sequencing libraries with poly A selection to enrich for protein coding transcripts. Then we carried our paired-end 2x150 Illumina sequencing. 

Our goal now is to quantify and compare gene expression between these two conditions by mapping and counting the numbers of sequencing reads that map to each gene in the reference transcriptome. 




