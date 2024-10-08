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

The output including read counts and measurement variance are contained in the `output` folder. Download this folder to you local computer, place it within a folder `DEG_analysis`, and we will continue the analyses within this folder locally using R. 

Before we leave the OSCER let's prepare a file that contains the transcript ID and the gene symbols. We will need this later. To do this we can use [grep](https://www.gnu.org/software/grep/manual/grep.html) to pull the headers from the transcriptome file and then we can pipe   `|` this to a [sed](https://www.gnu.org/software/sed/manual/sed.html) command that will parse the headers and write the ID and gene symbols to a new file. 

```
grep  '^>' GRCh38_latest_rna_with_ORF.fna | sed -E 's/>([^ ]+) .* \(([^)]+)\).*/\1 \2/' > TTC_headers.txt
```
Download `TTC_headers.txt` to the `DEG_analysis` folder on your local computer as well. 

## Differential gene expression analysis

To compare transcript expression levels between our treatment groups we will use [sleuth](https://pachterlab.github.io/sleuth/) a package designed specifically to work with kallisto output. We will also use the [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano) package for plotting. To install these, open up rstudio on your computer and run the following commands: 

```
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install("devtools")
BiocManager::install("pachterlab/sleuth")
```

In a text editor create tab-separated list of treatment conditions for our samples: 

```
Run_s	treat
CB1	control
CB2	control
CB3	control
CBT1	TTC
CBT2	TTC
CBT3	TTC
```
Save this as a text file `ExpTable_TTC.txt` in the `DEG_analysis` folder. Now that we have all of the files we need, set your working directory in R to `DEG_analysis`.

Now in R load the the packages we will be using: 

```
library(sleuth)
library(tidyverse)  
library(EnhancedVolcano)
library(pheatmap)
```
If you get a warning indicating the package is not installed, fo ahead and install through Rstudio's package manager (in the dropdown menu Tools > Install packages). 

Our next step is to setup the MetaData for our samples. 

```
#read in sample tables - be sure to set correct path 

metadata <- read.table(file = "ExpTable_TTC.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)

#this command sets up paths to the kallisto output that we will process in the following steps

metadata <- dplyr::mutate(metadata,
                          path = file.path('output', Run_s, 'abundance.h5'))
metadata <- dplyr::rename(metadata, sample = Run_s)

#let's check the metadata

> metadata
  sample   treat                     path
1    CB1 control  output/CB1/abundance.h5
2    CB2 control  output/CB2/abundance.h5
3    CB3 control  output/CB3/abundance.h5
4   CBT1     TTC output/CBT1/abundance.h5
5   CBT2     TTC output/CBT2/abundance.h5
6   CBT3     TTC output/CBT3/abundance.h5

#Read in headers for the transcripts that we aligned to with kallisto
#These will be mapped in the sleuth_prep command below

ttn<-read_delim("TTC_headers.txt", col_names = FALSE)

colnames(ttn)<-c("target_id","gene")
```
Now we are ready to process our data with the `sleuth_prep` command. This will aggregate the data from each samples individual result folder. 

```
so <- sleuth_prep(metadata, full_model = ~treat, target_mapping = ttn, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, aggregation_column = "gene")
```
`full_model` specifies the model we want to analyze, here we are interested in the effect of treatment. `target_mapping` will map our gene symbols to the transcript IDs. The boorstrap commands canculate stattistics for supplemental analyses we will explore with sleuth_live. `aggregation_column` sets up the dataset so we can aggregate transcripts of the same gene for our analyses. 

Nest, fit the model and calculate the test statisitics
```
#fit model specified above
so <- sleuth_fit(so)

#print the model
models(so)

# [  full  ]
# formula:  ~treat 
# data modeled:  obs_counts 
# transform sync'ed:  TRUE 
# coefficients:
# 	(Intercept)
#  	treatTTC

#calculate the Wald test statistic for 'beta' coefficient on every transcript 
so <- sleuth_wt(so, 'treatTTC')
```
To run the Wald test we need to specify which grouping to test. You can find the here `'treatTTC'` specifies the grouping by treatment. I simply took this from the model coefficients above. 

The results are now part of the sleuth object `so`. To recover the test results, we can run sleuth_results: 
```
#extract the wald test results for each transcript 
transcripts_all <- sleuth_results(so, 'treatTTC', show_all = FALSE, pval_aggregate = FALSE)
```
`pval_aggregate = FALSE` returns the transcript level comparisons, `show_all = FALSE` drops observations with missing data. 

Let's take a look: 
```
> head(transcripts_all, 10)
        target_id   gene          pval          qval         b      se_b mean_obs   var_obs
1  NM_001384818.1 PRRC2B 3.257746e-112 2.333751e-107  9.180827 0.4078420 3.897266 25.370898
2     NM_201378.4   PLEC  7.133144e-45  2.554985e-40  5.454296 0.3880557 4.709731  9.105510
3     NM_005560.6  LAMA5  2.119075e-38  5.060139e-34  8.168898 0.6304172 3.391302 20.046290
4     HOFI_TTC39B   <NA>  3.125429e-37  5.597409e-33  8.024288 0.6293657 3.948250 19.792079
5     NM_130440.4  PTPRF  4.094137e-34  5.865834e-30  1.309210 0.1075102 7.793013  0.526805
6     NM_014680.5  BLTP2  1.705087e-33  2.035788e-29  8.340445 0.6915440 3.746274 20.962528
7     NM_005334.3  HCFC1  9.665174e-33  9.891201e-29  7.978311 0.6694962 3.296008 19.105366
8  XM_005250983.3   PLEC  1.924058e-31  1.722922e-27  2.685866 0.2302500 6.928724  2.227780
9  XM_047421892.1   PLEC  5.797167e-31  4.614352e-27 -3.599663 0.3111001 5.099954  3.946645
10    NM_002184.4  IL6ST  8.606497e-29  6.165436e-25  7.793514 0.6999980 3.329170 18.315942
      tech_var     sigma_sq smooth_sigma_sq final_sigma_sq
1  0.002082568  0.103697515     0.247420124     0.24742012
2  0.111467572  0.114413321     0.070724895     0.11441332
3  0.014533950  0.019241107     0.581604909     0.58160491
4  0.048170876  0.545980961     0.226060846     0.54598096
5  0.011962814  0.003781762     0.005374840     0.00537484
6  0.392200834 -0.275173176     0.325148830     0.32514883
7  0.014576249 -0.002908329     0.657761394     0.65776139
8  0.048359599  0.031162979     0.008232621     0.03116298
9  0.104605526 -0.030388724     0.040569419     0.04056942
10 0.103492366  0.014362938     0.631503481     0.63150348
```
- `qval` is a false discovery rate adjusted p-value, using Benjamini-Hochberg method
- `b` effect size an estimator of the fold change between conditions

This is a huge table (>70,000 genes), but we are only really interested in the transcripts that are significantly differentially expressed. We can filter by qval using dplyr: 

```
#filtered by significance 
transcripts_sig <- dplyr::filter(transcripts_all, qval <= 0.05)
```
This yields a smaller list of ~1000 genes. The table is sorted by qval so we can select the 50 most significant transcripts by taking the first 50 entries. 

```
transcripts_50 <- dplyr::filter(transcripts_all, qval <= 0.05) %>%
  head(50)
```
Above, you may have noticed that there are multiple significant transcripts for the same gene (PLEC). Sleuth offers a way aggregate these transcripts and calculate gene-level tests. To do this we will set `pval_aggregate = TRUE`

```
genes_all <- sleuth_results(so, 'treatTTC', show_all = FALSE, pval_aggregate = TRUE)

> head(genes_all, 10)
   target_id num_aggregated_transcripts sum_mean_obs_counts          pval          qval
1       PLEC                         31           104.05922 3.247615e-208 7.959903e-204
2     PRRC2B                          4            24.04712 4.374874e-128 5.361408e-124
3      HCFC1                          9            25.87346  1.019680e-51  8.330782e-48
4     PTPN13                         17            72.39147  5.179964e-37  3.174023e-33
5      LAMA5                          5            18.15877  4.211684e-36  2.064568e-32
6      BLTP2                          7            28.34303  1.676552e-32  6.848713e-29
7      PTPRF                          7            27.65212  2.520152e-30  8.824131e-27
8    GOLGA8A                         15            74.21208  3.195410e-30  9.789938e-27
9     TRIP12                         35           108.74276  9.269410e-30  2.524369e-26
10      UBR4                         30           102.90138  4.706700e-27  1.153612e-23
```

### Volcano plots

A useful way to visualize the relative expression of all of the transcripts between two conditions is with a volcano plot that plots each transcript on an x-axis of fold-change in expression and y-axis of p-value. The most differential expressed genes will be found at the upper right and left cornerd of the plot. To do this we will use the EnhancedVolcano package that offers many ways to customize these plots.   

```
#extract the gene symbols, qval, and b values from the Wlad test results
forVolacano<-data.frame(transcripts_all$gene, transcripts_all$qval, transcripts_all$b)

#rename the columns of the dataframe
colnames(forVolacano)<-c("gene","qval","b")

#plot
EnhancedVolcano(forVolacano,
                lab = forVolacano$gene,
                x = 'b',
                y = 'qval',
                xlab = "\u03B2",
                labSize = 3,
                legendPosition = "none")
```

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image4.png)

### Heatmap of transcripts

Another common way to visualize differential gene expression among treatment groups and samples is with a heatmap. To plot the heat map we will need to exact the counts for each transcript and sample. This data is contained within the sleuth object `so` and we can export it with the `kallisto_table`
```
k_table <- kallisto_table(so, normalized = TRUE)
```
The `normalized` option will return values that have normalized for variation is sequencing depth and composition across the samples. Sleuth uses the [DESeq2 method](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html) for count normalization. We will plot the transcripts per million reads (tpm) which accounts for sequencing depth and gene length for each transcript. 

For plotting we need apply a log10 transfomation and convert this dataset into a square matrix: 
```
k_DEG_select<-k_DEG %>%
  #apply log10 transformation to the tpm data
  mutate(log_tpm = log10(tpm+1)) %>%
  #select the specifc columns to plot
  dplyr::select(target_id, sample, log_tpm, gene) %>%
  #create "label" from the transcript id and gene symbol
  mutate(label = paste(target_id, gene))%>%
  #pivot data frame to a wide format
  pivot_wider(names_from = sample, values_from = log_tpm) %>%
  #drop the target_id and gene variables
  dplyr::select(!target_id & !gene) %>%
  #convert label to row name
  column_to_rownames("label") %>%
  #convert to matrix
  as.matrix(rownames.force = TRUE) 

#plot with pheatmap!
pheatmap(k_DEG_select, cexRow = 0.4, cexCol = 0.4, scale = "none")
```
![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/RNAseq_image5.png)




