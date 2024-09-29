# Annotation - Gene finding and homology search

Genomes and transcriptomes are of little use for biological inferences if we do not know what is contain within them. However, annotation is a major challenge because ultimately our understanding of the function of sequence elements remains incomplete. Even for protein coding genes, which we understand the best, the validation of genes and transcripts in well-studied species remains a work in progress. Therefore the annotation process and approaches vary among fields and taxonomic groups. 


## Genome annotation

The current state of the art for eukaryotic genome annotation include [GNOMON](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/gnomon/) a proprietary analysis pipeline run by the NCBI and [Breaker](https://github.com/Gaius-Augustus/BRAKER) from the University of Greifswald. 

![GNOMON](https://www.ncbi.nlm.nih.gov/core/assets/genome/images/Pipeline_sm_ncRNA_CAGE_80pct.png)

These are complex multi-component pipelines involving dozens of separate software packages and integrating diverse data sources. These pipelines are beyond the scope of our workshop. However, we can implement key elements of these pipelines. 

All approaches involve some form of gene finding/prediction and homology searches with existing databases of the proteins. Here we will work through a simplified pipeline to find and annotate genes in [chromosome 24](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009792885.1/) of the horned lark

![Horned lark](https://www.allaboutbirds.org/guide/assets/photo/308604021-1280px.jpg)

### gene finding

To find genes in chromosome 24 we will use the package [Augustus](https://github.com/Gaius-Augustus/Augustus/tree/master) an ab initio program that bases its predictions on sequence data alone. Augustus relies on a [Hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) trained with RNA-seq, protein databases, an other data to calculate a probability that a particular sequence belongs to a gene. The creators offer trained models for many [common species](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/ABOUT.md) and there are also options to [train your own model](https://bioinf.uni-greifswald.de/webaugustus/trainingtutorial). Here we will use the existing model for chicken. 

```
ml AUGUSTUS/3.4.0-foss-2020b

augustus --species=chicken --protein=on /home/mbtoomey/BIOL7263_Genomics/Example_data/HoLa_scaffold_123.fasta > /scratch/mbtoomey/HoLa_example/HoLa_scaffold_123.gff

getAnnoFasta.pl /scratch/mbtoomey/HoLa_example/HoLa_scaffold_123.gff
```

The first command runs the gene finding with the chicken model and I have selected `--protein=on` so that augustus will out put predicted AA sequences. The `getAnnoFasta.pl` script writes these AA sequences to a new file that we will use next for homology search. 

This process yields two files a [.gff](https://useast.ensembl.org/info/website/upload/gff.html) that maps the predicted genes to the genome and a `.aa` file that contains the AA sequences of the predicted genes. 

* [augustus.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/augustus.sh)
* [augustus.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/augustus.sbatch)

### annotation by homology 

Now that we have putative protein coding genes in our genome we need to figure out what these might be. To do this we will compare them to known or annotated proteins in other bird species. We will use [Diamond](https://github.com/bbuchfink/diamond) to run efficient blast searches of databases. 

We are going to build our own database for this search. To do this lets go to the [Uniprotkb](https://www.uniprot.org/help/uniprotkb) database and download all of the protein sequences for chicken and zebra finch, two well annotated species. 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/annotation_image1.png)

I then simply unzipped these files and concatenated them into a single fasta:

```
cat ZeFi_proteins.fasta chicken_proteins.fasts > bird_proteins.fasta
```
Now we can convert out fast file of proteins into a blast database with diamond: 
```
diamond makedb --in /scratch/mbtoomey/HoLa_example/bird_proteins.fasta -d bird_proteins
```
* [diamond_mkdb.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb.sh)
* [diamond_mkdb.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb.sbatch)

Now we have a database file called `bird_proteins.dmnd` that we can search with our protein predictions from augustus. Let's do it!

```
diamond blastp --threads 20 --outfmt 6 -k 1 -d /scratch/mbtoomey/HoLa_example/bird_proteins.dmnd -q /scratch/mbtoomey/HoLa_example/HoLa_scaffold_123.aa -o HoLa_blastp.tsv
```

* `-d` sets the path to our database. 
* `-k` limits the output to the top hit only. 
* `-q` sets the path to our query sequences from the genome
* `-o` sets the name of the output file
* `--threads` sets the number of cpu cores to use for the anlysis. Remember to match this in the .sbatch file. 
* `--outfmt 6` sets the output to a table. By default, there are 12 fields included: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore. However this can be customized by specifying: 

     - `qseqid`
    Query Seq - id
    
    - `qlen`
    Query sequence length
    
    - `sseqid`
    Subject Seq - id
    
    - `sallseqid`
    All subject Seq - id(s), separated by a ’;’
    
    - `slen`
    Subject sequence length
    
    - `qstart`
    Start of alignment in query*
    
    - `qend`
    End of alignment in query*
    
    - `sstart`
    Start of alignment in subject*
    
    - `send`
    End of alignment in subject*
    
    - `qseq`
    Aligned part of query sequence*
    
    - `qseq_translated`
    Aligned part of query sequence (translated)* *Supported since v2.0.7.*
    
    - `full_qseq`
    Full query sequence
    
    - `full_qseq_mate`
    Query sequence of the mate (requires two files for `--query`) *Supported     since v2.0.7.*
    
    - `sseq`
    Aligned part of subject sequence*
    
    - `full_sseq`
    Full subject sequence
    
    - `evalue`
    Expect value
    
    - `bitscore`
    Bit score
    
    - `score`
    Raw score
    
    - `length`
    Alignment length*
    
    - `pident`
    Percentage of identical matches*
    
    - `nident`
    Number of identical matches*
    
    - `mismatch`
    Number of mismatches*
    
    - `positive`
    Number of positive - scoring matches*
    
    - `gapopen`
    Number of gap openings*
    
    - `gaps`
    Total number of gaps*
    
    - `ppos`
    Percentage of positive - scoring matches*
    
    - `qframe`
    Query frame
    
    - `stitle`
    Subject Title
    
    - `salltitles`
    All Subject Title(s), separated by a ’\<\>’

    - `qtitle`
    Query title
    
    - `qqual`
    Query quality values for the aligned part of the query*
    
    - `full_qqual`
    Query quality values
    

* [diamond_blastp.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastp.sh)
* [diamond_blastp.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastp.sbatch)

Now lets uses these blast hits to add features to genes found by augustus. Here we will edit the .gff file with [AGAT](https://agat.readthedocs.io/en/latest/index.html#) a convient tool for editing .gff files. There we some complicated dependencies so I set this up in a separate environment. You will need to activate it:
```
mamba activate /home/mbtoomey/.conda/envs/agat
```

```
agat_sp_manage_functional_annotation.pl -f HoLa_scaffold_123.gff -b HoLa_blastp.tsv --db bird_proteins.fasta --output HoLa_annotated.gff
```

This script takes our blastp output then searches the original protein database `bird_proteins.fasta` and pulls gene names from the corresponding headers and adds them to the annotation and makes a new annotation `HoLa_annotated.gff`

* [agat.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/agat.sh)
* [agat.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/agat.sbatch)

Now let's download our annotated GFF and genome file and fire up IGV. How did we do? Here is my favorite gene BCO2! :tada:

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/annotation_image1.png)




