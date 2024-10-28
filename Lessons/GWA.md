# Genome Wide Association

Genome-wide association (GWA) analyses aim to characterize the genetic differences among populations and identify genotype to phenotype associations. These studies involve the resequencing of individuals, mapping of sequence reads to a high-quality references genome, the calling of genetic variants (usually SNPs), and statistical approaches to test and characterize the distribution of those variants among samples.

For this example, we will replicate analyses from a [Aguillon et al. (2021)](https://royalsocietypublishing.org/doi/10.1098/rspb.2020.1805) of the northern flicker (*Colaptes auratus*). This species occurs as two subspecies inhabiting eastern North America (yellow-shafted) and western North America (red-shafted) with a zone of hybridization in the great plains. 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/GWAS_image1.png)

The aim of this study was to assess genomic divergence between the subspecies and analyze hybrid individuals to identify loci associated with plumage color phenotype. 

This is an exceptionally well documented study and I encourage you to check out the associated [Github repository](https://github.com/stepfanie-aguillon/flicker-WGS-ProcB2021/tree/master?tab=readme-ov-file). This is an excellent example of how you should aspire to organize your work!  

## Genome considerations

A successful GWA requires a high quality reference genome. Elsewhere we have covered [assembly](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_1.md) and [genome evaluation](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/Genome_Eval.md). The authors of this paper also scaffold their assembly with [Satsuma](https://github.com/bioinfologics/satsuma2) to generate a chromosomal level assembly and masked the repetative seqences with [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker). This is an important step because the repetitive sequences and associated assembly and alignment problems can generate spurious associations. 

## Alignment 

Once the reference is complete, the next step is to align the reads from each individual sample and output a bam file. This process is covered in detail in [Genome Adventure Chapter 2](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_6.md). 

## Variant calling

The next step is to call the variants (SNPs) from the bam files, this is covered in [Genome Adventure Chapter 2](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_14.md). In the flicker paper they used [GATK](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/) for their variant calls. 

## Set up for our exercise

For our exercise, we will start with the variant call file from Aguillon et al. (2021), but I have trimmed this down to just chromosome 3 to speed things up. 

For this exercise let's setup a folder in your scratch: 
```
cd /scratch/[your id]

mkdir flicker_gwa

cd flicker_gwa
```
Now make a symbolic link to the vcf file, genome, sample info:
```
ln -s /scratch/mbtoomey/flicker_gwa/flicker_chr3.recode.rename.vcf flicker_chr3.recode.rename.vcf
ln -s /scratch/mbtoomey/flicker_gwa/sample_info.txt sample_info.txt
ln -s /scratch/mbtoomey/flicker_gwa/NOFL_Chr3.fasta NOFL_Chr3.fasta
ln -s /scratch/mbtoomey/flicker_gwa/NOFL_Chr3.gff NOFL_Chr3.gff
```
To access the programs we will need for the next few steps activate the GWA environment:
```
mamba activate /home/mbtoomey/.conda/envs/GWA
```

## Genetic distance calculation between populations 

We can start by assessing the level of genetic differentiation between allopatric populations of the eastern, yellow-shafted and western, red-shafted flickers by calculating the fixation index (*F<sub>ST</sub>*). This a classic measure of differentation that compares the number of SNPS within and between populations. We can calculate this for the entire genome, but for our purposes we would like to scan the genome for regions with high *F<sub>ST</sub>* to do this we will define a window of 10 kb and then slide this along the chromosome in 2 kb steps. 

We will use vcftools to calculate the Weir and Cockerham [(1984)](https://www.jstor.org/stable/2408641) *F<sub>ST</sub>*. First, set up two text files listing the sample in each population. You can do this locally and upload, or try [nano](https://www.nano-editor.org/) in the command line.  

RSFL_indivs.txt:
```
samplenofl1toNOFL_sorted
samplenofl4toNOFL_sorted
samplenofl5toNOFL_sorted
samplenofl6toNOFL_sorted
samplenofl8toNOFL_sorted
samplenofl9toNOFL_sorted
samplenofl10toNOFL_sorted
samplenofl14toNOFL_sorted
samplenofl21toNOFL_sorted
samplenofl22toNOFL_sorted
```
YSFL_indivs.txt
```
samplenofl3toNOFL_sorted
samplenofl7toNOFL_sorted
samplenofl11toNOFL_sorted
samplenofl12toNOFL_sorted
samplenofl13toNOFL_sorted
samplenofl15toNOFL_sorted
samplenofl17toNOFL_sorted
samplenofl19toNOFL_sorted
samplenofl23toNOFL_sorted
samplenofl24toNOFL_sorted
```
Once these files are in place, we can submit a job to run the following: 

```
vcftools --vcf flicker_chr3.recode.rename.vcf --out RSFL_YSFL_chr3_10kb --weir-fst-pop RSFL_indivs.txt --weir-fst-pop YSFL_indivs.txt --fst-window-size 10000 --fst-window-step 2000`
```

`--fst-window-size` specifies the number of base pairs that will be used in the calculation and `--fst-window-step` specifies size of the step in base pairs as the analysis moves along the chromosomes. 

Here are my .sh and .sbatch files
* [fst.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/fst.sh)
* [fst.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/fst.sbatch)

The result is a the `RSFL_YSFL_10kb.windowed.weir.fst` file with F~ST~ values for each step: 

```
CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
Chr3    1               10000   71              -0.0178902      -0.0100243
Chr3    2001            12000   35              -0.0222068      -0.0145429
Chr3    4001            14000   35              -0.0222068      -0.0145429
Chr3    6001            16000   35              -0.0222068      -0.0145429
Chr3    8001            18000   37              -0.0248233      -0.0176512
Chr3    10001           20000   7               -0.070752       -0.067913
Chr3    12001           22000   7               -0.070752       -0.067913
Chr3    14001           24000   26              -0.0258256      -0.0257457
Chr3    16001           26000   26              -0.0258256      -0.0257457
```

Hang on to this file for now. We will come back to it after we complete the GWA analysis.

## Genome wide associations for feather color in hybrid individuals

To identify variants that are specifically associated with plumage coloration we will examine hybrid individuals that are the product of crosses between the red and yellow-shafted subspecies.

### Phasing

To prepare the variant data for GWA we will use [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) to:

* phase the genotypes (i.e. infer haplotypes) of the individuals in the samples.
* infer sporadic missing genotype data.

```
beagle gt=flicker_chr3.recode.rename.vcf nthreads=20 out=flicker_chr3.phased
```

* [beagle.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/beagle.sh)
* [beagle.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/beagle.sbatch)

### Filter the vcf to select hybrid individuals 

We want to focus our GWA on hybird inidviduals so referencing the sample data table for the ids, we will drop the allopatric individuals from the phased vcf file. 
```
vcftools --gzvcf flicker_chr3.phased.vcf.gz --remove-indv \
samplenofl1toNOFL_sorted --remove-indv\
samplenofl2toNOFL_sorted --remove-indv\
samplenofl3toNOFL_sorted --remove-indv\
samplenofl4toNOFL_sorted --remove-indv\
samplenofl5toNOFL_sorted --remove-indv\
samplenofl6toNOFL_sorted --remove-indv\
samplenofl7toNOFL_sorted --remove-indv\
samplenofl8toNOFL_sorted --remove-indv\
samplenofl9toNOFL_sorted --remove-indv\
samplenofl10toNOFL_sorted --remove-indv\
samplenofl11toNOFL_sorted --remove-indv\
samplenofl12toNOFL_sorted --remove-indv\
samplenofl13toNOFL_sorted --remove-indv\
samplenofl14toNOFL_sorted --remove-indv\
samplenofl15toNOFL_sorted --remove-indv\
samplenofl16toNOFL_sorted --remove-indv\
samplenofl17toNOFL_sorted --remove-indv\
samplenofl18toNOFL_sorted --remove-indv\
samplenofl19toNOFL_sorted --remove-indv\
samplenofl20toNOFL_sorted --remove-indv\
samplenofl21toNOFL_sorted --remove-indv\
samplenofl22toNOFL_sorted --remove-indv\
samplenofl23toNOFL_sorted --remove-indv\
samplenofl24toNOFL_sorted --recode --out flicker_chr3.phased.HZ
```

* filter_HZ.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/filter_HZ.sh)
* [filter_HZ.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/filter_HZ.sbatch)

### Prepare plink format files and add phenotype data 

Next we need to reformat the files with [Plink](https://zzz.bwh.harvard.edu/plink/) into .bed, .fam, and .bim formats for further analyses. 

`vcftools --gzvcf flicker_chr3.phased.HZ.recode.vcf --plink --out flicker_chr3_outputPlinkformat`

then

`plink --file flicker_chr3_outputPlinkformat --make-bed --out flicker_chr3_HZ_bed`

* [plink.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/plink.sh)
* [plink.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/plink.sbatch)

Now we should have three files: 

- `flicker_chr3_HZ_bed.bed`
- `flicker_chr3_HZ_bed.bim` 
- `flicker_chr3_HZ_bed.nosex`
- `flicker_chr3_HZ_bed.fam`


We now need to edit the `flicker_chr3_HZ_bed.fam` to add columns for the phenotype measures from `sample_info.txt`. To do this I download these two files to my local PC and ran the following in R: 

```
require(tidyverse)
fam<-read_table("flicker_chr3_HZ_bed.fam", col_names = FALSE)
fam<-fam[,1:5]
pheno<-read_table("sample_info.txt")

fam2<- pheno %>%
  select(ID, crown, ear, throat, nuchal, shaft) %>%
  left_join(fam, ., join_by(X1 == ID))

write_delim(fam2, "flicker_chr3_HZ_bed.fam", delim = " ", col_names = FALSE)
```

Now columns 6-10 of the .fam file will contain the phenotype categorizations. Later in gemma we will refer to these columns starting with column 6 = 1: 
```
# N1 = crown
# N2 = ear coverts
# N3 = throat
# N4 = nuchal patch
# N5 = wings and tail
```
Now upload the edited .fam file back to your `flicker_gwa` folder on OSCER. 

### Generate the relatedness matrix

The shared ancestry of the individuals in our samples will confound genome-type to phenotype relationships and should be carefully considered at the sampling stage. To model associations we need to provide measures of relatedness among our samples. Here we will infer relatedness from the variant dataset. 

For the GWA we will use [gemma](https://github.com/genetics-statistics/GEMMA) there are some conflict in the dependencies with the other pograms so we will run this in a separate environment:

```
mamba activate /home/mbtoomey/.conda/envs/gemma
```
Now you are set to run

```
gemma -bfile flicker_chr3_HZ_bed -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o GEMMA_HZ
```

* `-gk` specifies type of relatedness matrix to generate (1: centered matrix; 2: standardized matrix.)
* `-miss` missingness threshold - since we imputed with beagle this is not a concern
* `-maf` specify minor allele frequency threshold - we want to include all alleles regardless of frequency
* `-r2`  specify r-squared threshold 
* `-hwe` specify HWE test p value threshold - filters based on Hardy-Weinberg equilibrium - we do not want to filter and so set this to 0 

* [relate.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/relate.sh)
* [relate.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/relate.sbatch)

This command will generate a new folder `output` containing the `GEMMA_HZ.cXX.txt` relatedness matrix. 

### Test for associations 

We are now ready to test for associations between specific variants and the plumage color phenotype. To do this we fit a univariate linear model with our SNPs as a predictor of plumage color and calculate. We will calculate an effect size for each SNP, test against the null of no association, and correct p-values for multiple comparisons. 

```
gemma -bfile flicker_chr3_HZ_bed -k output/GEMMA_HZ.cXX.txt -lmm 4 -n 5 -o GWAS_HZ_lmm_shaft
```
* `-lmm 4` performs three tests: Wald test (-lmm 1), likelihood ratio test (-lmm 2), score test (-lmm 3)
* `-n 5` specifies the column containing feather color phenotype classifications.

* [gemma_lm.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/gemma_lm.sh)
* [gemma_lm.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/gemma_lm.sbatch)

The results are contained in the `GWAS_HZ_lmm_shaft.assoc.txt` file in the `output` folder: 

```
(gemma) [mbtoomey@schooner1 flicker_gwa]$ head output/GWAS_HZ_lmm_shaft.assoc.txt
chr     rs      ps      n_miss  allele1 allele0 af      beta    se      logl_H1 l_remle l_mle   p_wald  p_lrt   p_score
0       Chr3:434        434     0       A       G       0.219   4.203605e-01    2.879551e-01    -8.182207e+01   1.000000e+05    1.000000e+05 1.511376e-01    1.403847e-01    1.516829e-01
0       Chr3:462        462     0       T       C       0.281   4.237189e-01    2.779468e-01    -8.172606e+01   1.000000e+05    1.000000e+05 1.342400e-01    1.240241e-01    1.355205e-01
0       Chr3:534        534     0       C       A       0.198   3.168226e-01    3.372210e-01    -8.245277e+01   1.000000e+05    1.000000e+05 3.523763e-01    3.394966e-01    3.467556e-01
0       Chr3:554        554     0       T       G       0.146   1.002448e-01    4.517855e-01    -8.288326e+01   1.000000e+05    1.000000e+05 8.253850e-01    8.207368e-01    8.217876e-01
0       Chr3:565        565     0       G       A       0.083   -2.718738e-01   5.277283e-01    -8.277086e+01   1.000000e+05    1.000000e+05 6.088959e-01    5.992348e-01    6.022795e-01
0       Chr3:724        724     0       A       G       0.198   -3.053677e-01   3.289204e-01    -8.246340e+01   1.000000e+05    1.000000e+05 3.580495e-01    3.451883e-01    3.523309e-01
0       Chr3:736        736     0       G       A       0.219   -5.419099e-02   3.072384e-01    -8.289271e+01   1.000000e+05    1.000000e+05 8.607690e-01    8.570395e-01    8.578545e-01
0       Chr3:760        760     0       C       T       0.365   -9.722079e-02   2.873761e-01    -8.284930e+01   1.000000e+05    1.000000e+05 7.366726e-01    7.298188e-01    7.315533e-01
0       Chr3:799        799     0       T       C       0.490   -1.667614e-01   2.477951e-01    -8.267379e+01   1.000000e+05    1.000000e+05 5.043264e-01    4.928568e-01    4.973467e-01
```

### Viewing results in IGV

We will use [IGV](https://igv.org) to view the results and explore the loci that are differentiated among populations and associated with the red/yellow phenotype in the hybrid birds. To do this download: 

- `NOFL_Chr3.fasta`
- `NOFL_Chr3.gff`
- `RSFL_YSFL_10kb.windowed.weir.fst`
- `GWAS_HZ_lmm_shaft.assoc.txt`

to a folder on your local PC. We need to modify the format of the FST and GWA results to be readable by IGV and we can use R to do this: 

```
require(tidyverse)

#Read in the FST result file
fst<-read_tsv("RSFL_YSFL_chr3_10kb.windowed.weir.fst")

#Subset the columns for display in IGV
fst.igv<-fst %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST)

#write the subsetted file to tab separated file with the extension .igv
write_tsv(fst.igv, "RSFL_YSFL_chr3_10kb.windowed.weir.fst.igv")

#Read in the GWA result file
gwa<-read_tsv("GWAS_HZ_lmm_shaft.assoc.txt")

#relabel the chromosome. Gemma had converted our label "Chr3" to "0"
gwa$chr<-"Chr3"

#Subset the columns for display in IGV
gwa.igv<-gwa %>%
  select(chr, ps, rs, p_wald)

#Rename the columns specific to the GWAS format in IGV
colnames(gwa.igv)<-c("CHR","BP","SNP","P")

#write the subsetted file to tab separated file with the extension .igv
write_tsv(gwa.igv, "GWAS_HZ_lmm_shaft.assoc.gwas")
```

Now open up IGV ang go to Genomes > Load from file and load the `NOFL_Chr3.fasta`. Next go to File > Load file and load each: 

- `NOFL_Chr3.gff`
- `RSFL_YSFL_chr3_10kb.windowed.weir.fst.igv`
- `GWAS_HZ_lmm_shaft.assoc.gwas`

You should now have a display that looks somthing like this: 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/GWAS_image2.png)

Now we are ready to zoom in and explore the loci associate with the differentiation peaks and significant SNPs. 



