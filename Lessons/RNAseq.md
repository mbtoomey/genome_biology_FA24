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

However, we do not know how TTC39B is enhancing carotenoid metabolism. Therefore, we wanted to know if and how *TTC39B* expression alters gene expression within cells. To do this we cultured [HEK293 cells](https://en.wikipedia.org/wiki/HEK_293_cells) and transfected these cells with expression constructs for *CYP2J19* and *BDH1L* or *CYP2J19*, *BDH1L*, and *TTC39B*

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image3.png)

After 48 hours of expression, we harvested these cells, extracted RNA, and prepare short-read sequencing libraries with poly A selection to enrich for protein coding transcripts. Then we carried our paired-end 2x150 Illumina sequencing. 

Our goal now is to quantify and compare gene expression between these two conditions by mapping and counting the numbers of sequencing reads that map to each gene in the reference transcriptome. 

### The Data

For this exercise let's set up a folder on your scratch 

```
cd /scratch/[your id]

mkdir RNAseq_Example

cd RNAseq_Example
```

Now let's link to the RNA-Seq reads in my folder. I completed read quality assessment and trimming as we have done in other exercise, so you will not need to do these steps. 

```
ln -s /scratch/mbtoomey/RNAseq_example/CB1_S88_R1_001.fastq.gz CB1_S88_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CB2_S89_R1_001.fastq.gz CB2_S89_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CB3_S90_R1_001.fastq.gz CB3_S90_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CB1_S88_R2_001.fastq.gz CB1_S88_R2_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CB2_S89_R2_001.fastq.gz CB2_S89_R2_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CB3_S90_R2_001.fastq.gz CB3_S90_R2_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT1_S91_R1_001.fastq.gz CBT1_S91_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT2_S92_R1_001.fastq.gz CBT2_S92_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT3_S93_R1_001.fastq.gz CBT3_S93_R1_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT1_S91_R2_001.fastq.gz CBT1_S91_R2_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT2_S92_R2_001.fastq.gz CBT2_S92_R2_001.fastq.gz
ln -s /scratch/mbtoomey/RNAseq_example/CBT3_S93_R2_001.fastq.gz CBT3_S93_R2_001.fastq.gz
```
Now I recommend running everything from within the `RNAseq_Example` folder. This will allow us to keep the paths simple.

### Reference transcriptome

Since these are human cells I download the human reference transcription [GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz ](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) and added the coding sequences of  *CYP2J19*, *BDH1L*, and *TTC39B* genes from our expression constructs. 

## Pseudoalignement and quantification with Kallisto 

To generate read counts we will use [Kallisto](https://pachterlab.github.io/kallisto/) a program that uses a process called [pseudoalignment](https://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html) to rapidly find reads that match the reference targets and generates counts. Pseudoalignment involves breaking the reference and reads into k-mers and matching these two sets. 

### Building the reference transcriptome index

The first step in a Kallisto analysis is building a k-mer index of the references transcriptome. This can take a little while, so no need for you to run this, you can access the index I have built

```
ml kallisto/0.46.1-foss-2019a

kallisto index -i Human_RNA_ref.idx GRCh38_latest_rna_with_ORF.fna
```
Here `index` is the commmand to build an index and `-i` sepcifies the file name for the index. 

You can create a symbolic link to my index for your analyses: 

```
ln -s /scratch/mbtoomey/RNAseq_example/Human_RNA_ref.idx Human_RNA_ref.idx
```
### Transcript quantification

Now we are ready to pseudoalign and count transcripts. First let's make a folder to output the results. 

```
cd /scratch/[your id]/RNAseq_Example

mkdir output
```

Now we can setup the Kallisto quantification: 

```
ml kallisto/0.46.1-foss-2019a

kallisto quant -i Human_RNA_ref.idx -t 20 -o ./output/CB1 -b 50 CB1_val_1.fq.gz CB1_val_2.fq.gz
```
Here `quant` is the command to run quantification, `-i` specifies the the index, `-t` sets the number of cpus (be sure to match this in the .sbatch), `-o` set the output path (separate folder for each sample), and `-b` is the number of bootstrap repetitions that are used to calculate measurement variance used in downstream analyses in sleuth. 

We need to repeat this exact same analysis for all of our samples, so we can submit this to OSCER as an array job. To so this we will rewrite our command with placeholders that refer the to positions in an arguments file. 
```
ml kallisto/0.46.1-foss-2019a

kallisto quant -i Human_RNA_ref.idx -t 20 -o $1 -b 50 $2 $3
```

Here `$1` will refer to the first column of arguments files, `$2` the second and `$3` the third. Here is what the arguments file looks like: 

```
./output/CB1 CB1_S88_R1_001.fastq.gz CB1_S88_R2_001.fastq.gz
./output/CB2 CB2_S89_R1_001.fastq.gz CB2_S89_R2_001.fastq.gz
./output/CB3 CB3_S90_R1_001.fastq.gz CB3_S90_R2_001.fastq.gz
./output/CBT1 CBT1_S91_R1_001.fastq.gz CBT1_S91_R2_001.fastq.gz
./output/CBT2 CBT2_S92_R1_001.fastq.gz CBT2_S92_R2_001.fastq.gz
./output/CBT3 CBT3_S93_R1_001.fastq.gz CBT3_S93_R2_001.fastq.gz
```

It is simply a space separated table. Now if we specify an array in our .sbatch, OSCER will set up separate analyses for each of line in the arguments file: 

```
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=kallisto_%J_stdout.txt
#SBATCH --error=kallisto_%J_stderr.txt
#SBATCH --job-name=kallisto
#SBATCH --array=1-6
# 

bash kallisto_quant.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" kallisto_quant.args)
```

`--array` sets the number of seperate analyses and `$(sed -n "${SLURM_ARRAY_TASK_ID}p" kallisto_quant.args)` directs SLURM to fill in the `$` portions of the .sh file with the values in the .args file. 

* [kallisto_quant.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/kallisto_quant.sh)
* [kallisto_quant.args](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/kallisto_quant.args)
* [kallisto_quant.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/kallisto_quant.sbatch)

If you take a look at the stderr files you will see details of the aligment process including the number of reads that pseudoaligned and estimate fragment lengths

```
The following have been reloaded with a version change:
  1) binutils/2.38 => binutils/2.31.1-GCCcore-8.2.0


[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 185,124
[index] number of k-mers: 142,407,651
[index] number of equivalence classes: 657,658
[quant] running in paired-end mode
[quant] will process pair 1: CBT1_S91_R1_001.fastq.gz
                             CBT1_S91_R2_001.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 45,331,634 reads, 33,134,805 reads pseudoaligned
[quant] estimated average fragment length: 222.854
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 1,517 rounds
[bstrp] number of EM bootstraps complete: 1
[bstrp] number of EM bootstraps complete: 2
[bstrp] number of EM bootstraps complete: 3
```

The output including read counts and measurement variance are contained in the output folder. Download this to you local computer and we will continue the analyses localy using R. 

## Differential gene expression analysis





