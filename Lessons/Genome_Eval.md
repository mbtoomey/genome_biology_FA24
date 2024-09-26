# Genome Assembly Evaluation

## Assembly statitistics

* **Genome size** - total # of bases in genome assembly
* **Contig/scaffold numbers** - total number of assembled elements
ideally = # of chromosomes
* **N50** - 50% of nucleotides in genome are contained in contigs that are greater than or equal to this length
* **L50** - minimum number of contigs needed to contain 50% of the genome

### Quast

[Quast - QUality ASsessment Tool](https://github.com/ablab/quast) is a convenient tool to assess the quality of your assembly and will report the statistics above and other details. 

Let's take a look at the pseudomonas genome assembly you will construct in Chapter 5 of the Genomics Adventure.

```
quast.py --output-dir /scratch/mbtoomey/pseud_quast_illumina/ /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/illumina/contigs.fasta
```
* [pseud_quast_illumina.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_illumina.sh)
* [pseud_quast_illumina.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_illumina.sbatch)

We simply provide our fasta file containing the assembled contigs and specify an output folder where the report documents will be saved. You can then view the report.txt in the terminal `cat /scratch/mbtoomey/pseud_quast_illumina/report.txt`

```
Assembly                    contigs
# contigs (>= 0 bp)         528
# contigs (>= 1000 bp)      117
# contigs (>= 5000 bp)      100
# contigs (>= 10000 bp)     89
# contigs (>= 25000 bp)     74
# contigs (>= 50000 bp)     51
Total length (>= 0 bp)      6692019
Total length (>= 1000 bp)   6618751
Total length (>= 5000 bp)   6580437
Total length (>= 10000 bp)  6497588
Total length (>= 25000 bp)  6248274
Total length (>= 50000 bp)  5443414
# contigs                   127
Largest contig              248605
Total length                6626613
GC (%)                      59.00
N50                         102161
N75                         64662
L50                         22
L75                         43
# N's per 100 kbp           0.00
```

or download the report folder to your local PC and view the html report files:

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_1.png)

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_2.png)

#### Assembly comparisons

Quast makes comparing different genome assemblies, of the same genome, easy. Simply provide as many contig fasta files as you would like to compare.

Here we can compare the assembly we will make with illumina short reads only to a hybrid assembly built with short-reads and long-reads: 

```
quast.py --output-dir /scratch/mbtoomey/pseud_quast_short_long/ /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/illumina/contigs.fasta /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/hybrid/contigs.fasta
```

* [pseud_quast_short_long.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_short_long.sh)
* [pseud_quast_short_long.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_short_long.sbatch)

Download the report folder and take a look at the report: 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_3.png)

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_4.png)

#### Reference comparison

Finally, is a high quality assembly is available you can compare your assemblies to it and quast will identify various misassemblies. 

Here I have downloaded the [Genome assembly Pseudomonas.strGM41_v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000282315.2/) from NCBI for a reference. 

```
quast.py --output-dir /scratch/mbtoomey/pseud_comp_quast/ /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/illumina/contigs.fasta /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/hybrid/contigs.fasta -r /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/pseud_ref.fna.gz
```

The `-r` option here specifies which file is the reference. 

* [pseud_quast_comp.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_short_long.sh)
* [pseud_quast_comp.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud_quast_short_long.sbatch)

Download the report folder and take a look at the report: 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_5.png)

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_6.png)

## Visualize the assembly graph

A graphical representation of your assembly can be helpful to: 

* Understand the connections among your contigs
* Diagnose problematic regions of the assembly
* Visualize location of genes and other regions of interest

[Bandage](https://rrwick.github.io/Bandage/) is a popular tool to visualize the assembly graph file and can be installed an run locally on your PC. 

Bandage requires a graph file. This is a file produced by your assembler. For example, SPAdes produces `.fastg` graph files along with the contig fasta. 

Here are the assembly graph files for pseudomonas genome assemblies you will construct in Chapter 5 of the Genomics Adventure. Let's download this and load them into Bandage to explore it's functions.

* [illumina assembly graph](https://drive.google.com/file/d/1RSAtTPTfF1nM0hos3u8elQrETdM0jNjF/view?usp=sharing)
* [hybrid assembly graph](https://drive.google.com/file/d/1Uxw-jzetJhYPOH-xHHzJqKqMYKYs-wh6/view?usp=sharing)

Let's try and find the [Universal Stress Protein](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/USP.fasta) in our assembly. 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_7.png)

## K-mer spectrum analysis

K-mer spectrum analysis is a useful tool to: 
* Assess quality of your sequencing library
* Infer properties of the genome
  * Size
  * Ploidy 
  * Heterozygosity
  
To do this we will use [Jellyfish](https://www.genome.umd.edu/jellyfish.html#Release) for K-mer counting and visualize the spectra with [Genomescope](http://genomescope.org/genomescope2.0/).

Let's take a look at the illumina sequencing reads we used to construct the pseudomonas genome assembly in Chapter 5 of the Genomics Adventure.

First we will count the k-mers and export a histogram file with Jellyfish

```
ml Jellyfish/2.3.0-GCCcore-7.3.0

jellyfish count -C -m 31 -s 31000000000 -t 20 <(zcat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_1.fastq.gz) <(zcat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_2.fastq.gz) -o /scratch/mbtoomey/pseud_merqury/pseud.jf

jellyfish histo -t 20 /scratch/mbtoomey/pseud_merqury/pseud.jf > /scratch/mbtoomey/pseud_merqury/pseud.jf.histo
```
Here I am doing a few things differently that usual. Since Jellyfish is installed on OSCER for all users we can activate it with the module load (ml) command. `-C` specified "canonical k-mers", `-m` sets the k-mer length, `-s` sets the amount of memory the program will use, and `-t` set the number of cores used. This is rather computationally intensive so I set it to use 20 cores. 

The second command `histo` generates the file that we will upload to genomescope. 

* [jellyfish.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/jellyfish.sh)
* [jellyfish.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/jellyfish.sbatch)

Let's upload the resulting [pseud.jf.histo](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/pseud.jf.histo) to Genomescope: 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/genome_eval_8.png)

## BUSCO Analysis to assess genome completeness

Benchmarking Universal Single-Copy Orthologs, [BUSCO](https://busco.ezlab.org/) is an approach that searches your assembly for a curated set of protein coding genes that are nearly universally present in the genomes of a selected taxonomic divison. 

Let's evaluate the pseudomonas genome we assembled: 

```
ml BUSCO/5.2.2-foss-2020b

busco -i /home/mbtoomey/BIOL7263_Genomics/precomp_hybrid_assembly/hybrid/contigs.fasta -m genome --lineage_dataset pseudomonadales_odb10 -c 20 --out pseud_busco2
```
`-i` specifies the input file - a .fasta of you contigs. `-m` sets the mode we are using genome, but you can also run transcriptome or proteins. `-c` set the number of CPU cores that busco will use. `--out` specifies the folder where results will be output. :heavy_exclamation_mark: you cannot specify a path here. The output folder will be created wherever you run the busco commands from. Also busco will create a database folder wherever you are running the command from. 

`--lineage_dataset` selects the taxonomic grouping that busco will search within. `busco --list-datasets` will output a list of all available lineages. If you are unsure what lineage to you you can instead specify `--auto-lineage` and busco will search your assembly and identify the most approporiate lineage, however this will loner than a specific search. 

* [busco.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/busco.sh)
* [busco.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/busco.sbatch)

```
# BUSCO version is: 5.2.2 
# The lineage dataset is: pseudomonadales_odb10 (Creation date: 2024-01-08, number of genomes: 159, number of BUSCOs: 782)
# Summarized benchmarking in BUSCO notation for file /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/contigs.fasta
# BUSCO was run in mode: genome
# Gene predictor used: prodigal

	***** Results: *****

	C:99.8%[S:99.7%,D:0.1%],F:0.1%,M:0.1%,n:782	   
	781	Complete BUSCOs (C)			   
	780	Complete and single-copy BUSCOs (S)	   
	1	Complete and duplicated BUSCOs (D)	   
	1	Fragmented BUSCOs (F)			   
	0	Missing BUSCOs (M)			   
	782	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.3
	prodigal: 2.6.3
	
```
	
	All in all, our assembly is pretty good! :clap:


