# Annotation - Gene finding and homology search

Genomes and transcriptomes are of little use for biological inferences if we do not know what is contain within them. However, annotation is a major challenge because ultimately our understanding of the function of sequence elements remains incomplete. Even for protein coding genes, which we understand the best, the validation of genes and transcripts in well-studied species remains a work in progress. Therefore the annotation process and approaches vary among fields and taxonomic groups. 


## Genome annotation

The current state of the art for eukaryotic genome annotation include [GNOMON](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/gnomon/) a proprietary analysis pipeline run by the NCBI and [Breaker](https://github.com/Gaius-Augustus/BRAKER) from the University of Greifswald. 

![GNOMON](https://www.ncbi.nlm.nih.gov/core/assets/genome/images/Pipeline_sm_ncRNA_CAGE_80pct.png)

These are complex multi-component pipelines involving dozens of separate software packages and integrating diverse data sources. These pipelines are beyond the scope of our workshop. However, we can implement key elements of these pipelines. 

All approaches involve some form of gene finding/prediction and homology searches with existing databases of the proteins. Here we will work through a simplified pipeline to find and annotate genes in [chromosome 24](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009792885.1/) of the horned lark

![Horned lark](https://www.allaboutbirds.org/guide/assets/photo/308604021-1280px.jpg)

### Gene finding

For this exercise let's set up a folder on your scratch 

```
cd /scratch/[your id]

mkdir HoLa_annotation

cd HoLa_annotation
```

Now let's link to the genome assembly and protein database stored in my home directory. To do this you will create a [symbolic link](https://servicenow.iu.edu/kb?id=kb_article_view&sysparm_article=KB0023928)

```
 ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/HoLa_scaffold_123.fasta HoLa_scaffold_123.fasta
 
ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/bird_proteins.fasta bird_proteins.fasta
```
Now I recommend running everthing from within the HoLA_example folder. This will allow us to keep the paths simple.

To find genes in chromosome 24 we will use the package [Augustus](https://github.com/Gaius-Augustus/Augustus/tree/master) an ab initio program that bases its predictions on sequence data alone. Augustus relies on a [Hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) trained with RNA-seq, protein databases, an other data to calculate a probability that a particular sequence belongs to a gene. The creators offer trained models for many [common species](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/ABOUT.md) and there are also options to [train your own model](https://bioinf.uni-greifswald.de/webaugustus/trainingtutorial). Here we will use the existing model for chicken. 

```
ml AUGUSTUS/3.4.0-foss-2020b

augustus --species=chicken --protein=on HoLa_scaffold_123.fasta > HoLa_scaffold_123.gff

getAnnoFasta.pl HoLa_scaffold_123.gff
```

The first command runs the gene finding with the chicken model and I have selected `--protein=on` so that augustus will out put predicted AA sequences. The `getAnnoFasta.pl` script writes these AA sequences to a new file that we will use next for homology search. 

This process yields two files a [.gff](https://useast.ensembl.org/info/website/upload/gff.html) that maps the predicted genes to the genome and a `.aa` file that contains the AA sequences of the predicted genes. 

* [HoLa_augustus.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/HoLa_augustus.sh)
* [HoLa_augustus.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/HoLa_augustus.sbatch)

### Annotation by homology 

Now that we have putative protein coding genes in our genome we need to figure out what these might be. To do this we will compare them to known or annotated proteins in other bird species. We will use [Diamond](https://github.com/bbuchfink/diamond) to run efficient blast searches of databases. 

We are going to build our own database for this search. To do this lets go to the [Uniprotkb](https://www.uniprot.org/help/uniprotkb) database and download all of the protein sequences for chicken and zebra finch, two well annotated species. 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/annotation_image1.png)

I then simply unzipped these files and concatenated them into a single fasta:

```
cat ZeFi_proteins.fasta chicken_proteins.fasta > bird_proteins.fasta
```
You will access this file through the symbolic link we setup above. Now we can convert our fasta file of proteins into a blast database with diamond: 
```
diamond makedb --in bird_proteins.fasta -d bird_proteins
```
* [diamond_mkdb.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb.sh)
* [diamond_mkdb.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb.sbatch)

Now we have a database file called `bird_proteins.dmnd` that we can search with our protein predictions from augustus. Let's do it!

```
diamond blastp --threads 8 --outfmt 6 -k 1 -d bird_proteins.dmnd -q HoLa_scaffold_123.aa -o HoLa_blastp.tsv
```

* [diamond_blastp.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastp.sh)
* [diamond_blastp.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastp.sbatch)

I have set the following options for this analysis: 

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

Now let's use these blast hits to add features to genes found by augustus. Here we will edit the .gff file with [AGAT](https://agat.readthedocs.io/en/latest/index.html#) a convient tool for editing .gff files. There we some complicated dependencies so I set this up in a separate environment. You will need to activate it:
```
mamba activate /home/mbtoomey/.conda/envs/agat
```

```
agat_sp_manage_functional_annotation.pl -f HoLa_scaffold_123.gff -b HoLa_blastp.tsv --db bird_proteins.fasta --output HoLa_annotated.gff
```

This script takes our blastp output then searches the original protein database `bird_proteins.fasta` and pulls gene names from the corresponding headers and adds them to the annotation and makes a new annotation `HoLa_annotated.gff`

* [agat.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/agat.sh)
* [agat.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/agat.sbatch)

Now let's download our annotated GFF and genome file and fire up IGV. To get the genome file you will need now need to copy the file to your sc

How did we do? Here is my favorite gene BCO2! :tada:

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/annotation_image2.png)

## Transcriptome annotation

We can use a similar approach to annotate **de novo** assembled transcriptomes produced from RNA sequencing data. For this example, I used Illumina PE150 reads from an RNA-seq project investigating gene expression in [HEK293 cells](https://en.wikipedia.org/wiki/HEK_293_cells). First I trimmed the reads and then assembled with [SPAdes in rna mode](https://ablab.github.io/spades/rna.html): 

```
spades.py --rna -t 20 -m 60 -o spades_assembly -1 trimmed_reads_val_1.fq.gz -2 trimmed_reads_val_2.fq.gz
```
This assembly required about 2 hrs, so we will not run it here, Rather you can access the finished assembly by setting a symbolic link to the file in my home directory. 

First lets make a folder for this exercise:
```
cd /scratch/[your id]

mkdir HEK_annotation

cd HEK_annotation
```

As above, I recommend running all of the scripts from this folder. Now let's make the links to the assembly and the human protein reference database:

```
ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/transcripts.fasta transcripts.fasta
 
ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/human_proteins.faa human_proteins.fasta
```

Similar to above, I download the human proteins from [Uniprotkb](https://www.uniprot.org/help/uniprotkb). Now we need to transform this file into a database that diamond can use for our blast search: 

```
diamond makedb --in human_proteins.fasta -d human_proteins
```

* [diamond_mkdb_HEK.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb_HEK.sh)
* [diamond_mkdb_HEK.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_mkdb_HEK.sbatch)

Now we can run the blast search, however this time we will be searching RNA transcripts against a protein database. Therefore, we will use the `blastx` search option in diamond. 

```
diamond blastx --threads 8 --outfmt 6 qseqid sseqid length pident evalue stitle -k 1 -d human_proteins.dmnd -q transcripts.fasta -o HEK_blastx.tsv
```

* [diamond_blastx.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastx.sh)
* [diamond_blastx.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/diamond_blastx.sbatch)

Now we have our blast results in tabular form `HEK_blast.tsv`. Let's merge these with out transcriptome assembly and add annotation to the tranacript headers. To do this lets first strip out the information we are interested in with an [awk](https://en.wikipedia.org/wiki/AWK) command. 

```
awk '{split($0, a, " "); split(a[2], b, "|"); gene = ($0 ~ /GN=/) ? gensub(/.*GN=([^ ]+).*/, "\\1", "g") : "-"; desc_start = index($0, a[7]); desc_end = match(substr($0, desc_start), / OS=/); print a[1] "\t" gene "|" a[2] " " substr($0, desc_start, desc_end - 1)}' HEK_blastx.tsv > HEK_headers.txt
```
This will take the blast output:
```
NODE_304989_length_197_cov_4.048387_g297389_i0  tr|E9PJ07|E9PJ07_HUMAN  21      100     2.52e-09        tr|E9PJ07|E9PJ07_HUMAN Erythroid differentiation regulatory factor 1 OS=Homo sapiens OX=9606 GN=EDRF1 PE=1 SV=1
NODE_304994_length_197_cov_4.032258_g297394_i0  tr|H0Y5H2|H0Y5H2_HUMAN  55      98.2    2.73e-28        tr|H0Y5H2|H0Y5H2_HUMAN Cdk5 and Abl enzyme substrate 2 (Fragment) OS=Homo sapiens OX=9606 GN=CABLES2 PE=1 SV=1
NODE_304996_length_197_cov_4.032258_g297396_i0	tr|Q6ZNN6|Q6ZNN6_HUMAN	26	65.4	9.86e-05	tr|Q6ZNN6|Q6ZNN6_HUMAN cDNA FLJ27422 fis, clone WMC08087 OS=Homo sapiens OX=9606 PE=2 SV=1
```
and parse it to this: 
```
NODE_304989_length_197_cov_4.048387_g297389_i0	EDRF1|tr|E9PJ07|E9PJ07_HUMAN Erythroid differentiation regulatory factor 1
NODE_304994_length_197_cov_4.032258_g297394_i0	CABLES2|tr|H0Y5H2|H0Y5H2_HUMAN Cdk5 and Abl enzyme substrate 2 (Fragment)
NODE_304996_length_197_cov_4.032258_g297396_i0	-|Q6ZNN6_HUMAN cDNA FLJ27422 fis, clone WMC08087
```

for the sake of transparency I didn't code this awk command myself, but rather submited [this query ](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/query.txt) to [copilot](copilot.microsoft.com) and then iterated until I got what I wanted. 

Now that we have the header file we can merge this with our transcriptome using the replace function in [seqkit](https://bioinf.shenwei.me/seqkit/) a package with many useful functions to edit fasta files. 

```
seqkit replace -p "(.+)" -r '$1|{kv}' -k HEK_headers.txt transcripts.fasta > transcripts_annotated.fasta
```

- `-p "(.+)"` This is the pattern to match the entire header. (.+) captures the entire header as a group.
- `-r '$1|{kv}'` This is the replacement pattern. $1 refers to the entire matched header (captured by (.+)), and {kv} is a placeholder that will be replaced by the corresponding value from the key-value file.
- `-k HEK_headers.txt` This specifies the key-value file. The file HEK_headers.txt should contain pairs of original headers and their replacements.
- `transcripts.fasta` This is the input FASTA file whose headers you want to replace.

The result is a fasta with the matching gene details added. Let's compare the original transcriptome fasta to the new annotated one:

```
grep '^>' transcripts.fasta | tail

>NODE_304997_length_197_cov_4.024194_g297397_i0
>NODE_304998_length_197_cov_4.016129_g297398_i0
>NODE_304999_length_197_cov_4.016129_g297399_i0
>NODE_305000_length_197_cov_4.008065_g297400_i0
>NODE_305001_length_197_cov_4.008065_g297401_i0
>NODE_305002_length_197_cov_4.008065_g297402_i0
>NODE_305003_length_197_cov_4.008065_g297403_i0
>NODE_305004_length_197_cov_4.008065_g297404_i0
>NODE_305005_length_197_cov_2.451613_g297405_i0
>NODE_305006_length_166_cov_2643.644351_g297406_i0

grep '^>' transcripts_annotated.fasta | tail

>NODE_304997_length_197_cov_4.024194_g297397_i0|
>NODE_304998_length_197_cov_4.016129_g297398_i0|ZNF568|tr|K7EKY2|K7EKY2_HUMAN Zinc finger protein 568
>NODE_304999_length_197_cov_4.016129_g297399_i0|COA8|tr|A0A0U1RR29|A0A0U1RR29_HUMAN Cytochrome c oxidase assembly factor 8 (Fragment)
>NODE_305000_length_197_cov_4.008065_g297400_i0|
>NODE_305001_length_197_cov_4.008065_g297401_i0|
>NODE_305002_length_197_cov_4.008065_g297402_i0|-|tr|Q7KZ41|Q7KZ41_HUMAN RNA-directed DNA polymerase (Fragment)
>NODE_305003_length_197_cov_4.008065_g297403_i0|TNNI2|sp|P48788|TNNI2_HUMAN Troponin I, fast skeletal muscle
>NODE_305004_length_197_cov_4.008065_g297404_i0|-|tr|Q6ZVR1|Q6ZVR1_HUMAN cDNA FLJ42200 fis, clone THYMU2034647
>NODE_305005_length_197_cov_2.451613_g297405_i0|
>NODE_305006_length_166_cov_2643.644351_g297406_i0|YWHAE/FAM22A|tr|G9K388|G9K388_HUMAN 14-3-3 protein epsilon (Fragment)
```
Here I used [grep](https://www.gnu.org/software/grep/manual/grep.html) to return just the headers from the file then piped `|` these to `tail` to look at the last few entries. 

Note that many of the assembled transcripts had no hits in our protein database. This is not surprising, because the transcriptome contains many RNAs not represented in this dataset (i.e. non-coding RNAs)

## BLAST

Diamond only allows for searches of a protein database. To search a nucleotide database we can use the [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) suite of tools. This is the exact same tool that is used for [online blast searches](https://blast.ncbi.nlm.nih.gov/Blast.cgi). However, we we implement it on our system we will need to provide a search database. The blast+ suit includes a script to download the preconfigured databases from NCBI: 

```
update_blastdb.pl --decompress [DATABASE NAME]
```
You can find details about how to use the script [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/)

To check the databases that are available you can run: 
```
update_blastdb.pl --showall

18S_fungal_sequences
Betacoronavirus
28S_fungal_sequences
ITS_RefSeq_Fungi
16S_ribosomal_RNA
ITS_eukaryote_sequences
LSU_eukaryote_rRNA
LSU_prokaryote_rRNA
SSU_eukaryote_rRNA
env_nt
env_nr
human_genome
landmark
mito
mouse_genome
nr
nt_euk
nt
nt_others
nt_prok
nt_viruses
pataa
patnt
pdbaa
pdbnt
ref_euk_rep_genomes
ref_prok_rep_genomes
ref_viroids_rep_genomes
ref_viruses_rep_genomes
refseq_select_rna
refseq_select_prot
refseq_protein
refseq_rna
swissprot
tsa_nr
tsa_nt
taxdb
core_nt
```

`core_nt` is the default database used in online blast and is >160 Gb compressed, so it is not practical to download the entire database. However, `nt_viruses` and `nt_prok` (prokaryote) databases are subsets that more reasonable in size and may be worth a try. Unfortunately, the eukaryotic database (`nt_euk`) is still impractically large for our purposes, so we will use a different approach to make our database for this example. 

To create a database of RNA sequences to search, I downloaded all of the RNAs from the NCBI human refernce genome - [GCF_000001405.40_GRCh38.p14_rna.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz). The I converted this to database:

```
makeblastdb -in human_RNA.fna -parse_seqids -blastdb_version 5 -title "Human RNA" -dbtype nucl -out human_rna_db
```

You can access with database in my account at ```/home/mbtoomey/BIOL7263_Genomics/Example_data/blastdb/human_rna_db```

Now we can run a nucleotide to nucleotide `blastn` search of our **de novo** transcriptome.
```
blastn -db /home/mbtoomey/BIOL7263_Genomics/Example_data/blastdb/human_rna_db -query transcripts.fasta -outfmt "6 qseqid sseqid stitle" -num_threads 20 -num_alignments 1 > TTC_rna_blast.tsv
```
The options here are a bit different than the diamond search:
- `-db` points toward our blast database (search target)
- `-query` is our assembly
- `num_threads` sets the number of CPU cores
- `num_aligments` sets the number of blast hits to return for each quary
- `outfmt` Here we have set it to tablular "6" and sepcified the specific elements that we want in the table. By default, there are 12 fields included: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore. However this can be customized by specifying: 

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
    
Here are my .sh an .sbatch to run the blast search

* [blastn.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/blastn.sh)
* [blastn.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/blastn.sbatch)

Now download and compare the results of the diamond protein search `HEK_blastx.tsv` to the blastn nucleotide search `HEK_blastn.tsv`. You will notice they are similar, but not identical. The blastn search identified more of the **de novo** transcripts as you might expect since this database contains non-coding RNAs and non-coding protions of the transcript sequences. Thus, a nucleotide-based search my be better, but keep in mind that nucleotide sequences are mush less conserved than protein sequences. If you are comparing to a database from a distantly releated taxa, the blastx protein search may be more reliable. 

### Side quest - blasting raw reads

A situation may arise where you would like to search your raw sequencing reads for a specific gene or sequence. Basic BLAST is not well suited to this task. However, there is an alternative version called [Magic-BLAST](https://ncbi.github.io/magicblast/) that is designed to query a database with sequencing reads. 

Let's work through an example with the HEK cell experiment raw reads. In this experiment I expressed several carotenoid metabolizing genes and was interested in if and how the exprerssion of these heterologous genes affected gene expression in the cells. One of those genes was [CYP2J19](), let's see if we can find this transcript in the raw reads.  You can access the CYP2J19 sequence and the raw reads from my directory. Now set up a folder for this example and create symbolic links to the files in my account. Note, to speed things up I subsampled the raw reads with seqtk as we did in [Cahpter 2](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_3.md) of the genomics adventure. 

```
cd /scratch/[your id]/

mkdir raw_read_blast

cd raw_read_blast

ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/CYP2J19.fasta CYP2J19.fasta

ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/subsample_1.fq subsample_1.fq

ln -s /home/mbtoomey/BIOL7263_Genomics/Example_data/subsample_2.fq subsample_2.fq
```

First we will need to transform the `CYP2J19.fasta` into a searchable database. As we did above, we will use the `makeblastdb` function. I installed magicblast to a separate environment, so you will first need to activate the `magicbalst` environment with mamba and then you can proceed. 

```
mamba activate /home/mbtoomey/.conda/envs/magicblast

makeblastdb -in CYP2J19.fasta -dbtype nucl -parse_seqids -out CYP2J19 -title "CYP2J19"
```
Now you should have a very small database, consisting of just the CYP2J19 nucleotide sequence

```
magicblast -query subsample_1.fq -query_mate subsample_1.fq -db CYP2J19 -infmt fastq -outfmt sam -no_unaligned -out HEK_CYP2J19_blast.sam
```
The options here are similar to basic blast:
- `-db` points toward our blast database (search target)
- `-query` and `-querymate` specify the forward and reverse of paired reads
- `infmt` specifies the format that the raw reads are in
- `num_aligments` sets the number of blast hits to return for each quary
- `outfmt` you can select a "sam" file or a table. More details [here](https://ncbi.github.io/magicblast/doc/output.html)
- `no_unaligned` returns only reads that align to the database, otherwise it will return all of the reads and mark aligned and unaligned with the sam flag
- `-out` specified the output file

* [mblast.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/mblast.sh)
* [mblast.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/mblast.sbatch)

The result of our search is a sam file containing all of the reads that matched the database: 

```
(base) [mbtoomey@schooner3 raw_read_blast]$ cat HEK_CYP2J19_blast.sam
@HD     VN:1.0  GO:query
@SQ     SN:HOFI_CYP2J19 LN:1494
@PG     ID:magicblast   PN:magicblast   CL:magicblast -query subsample_1.fq -query_mate subsample_1.fq -db CYP2J19 -infmt fastq -outfmt sam -no_unaligned -out HEK_CYP2J19_blast.sam
LH00260:42:222HHLLT4:2:1136:6300:19728  113     HOFI_CYP2J19    241     60      151M    =       241     -151CAGTTTGGAAGTCTGACATTCGTGGTGGTCAACGGGTACCAGATGGTGAGAGAAGCTCTTGTCCACCAGGCTGAAATATTTGCTGACCGGCCAAATATTCCACTCCTCCAAGAAATATTTAGAGGCTTTGGGCTCATATCATCCAACGGGC     II9IIIIIIIIIII9IIIIIIII9IIIII9IIIIIIIIIIIIII9IIIIIII-IIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     NH:i:1  AS:i:151    NM:i:0
LH00260:42:222HHLLT4:2:1136:6300:19728  177     HOFI_CYP2J19    241     60      151M    =       241     -151CAGTTTGGAAGTCTGACATTCGTGGTGGTCAACGGGTACCAGATGGTGAGAGAAGCTCTTGTCCACCAGGCTGAAATATTTGCTGACCGGCCAAATATTCCACTCCTCCAAGAAATATTTAGAGGCTTTGGGCTCATATCATCCAACGGGC     II9IIIIIIIIIII9IIIIIIII9IIIII9IIIIIIIIIIIIII9IIIIIII-IIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     NH:i:1  AS:i:151    NM:i:0
LH00260:42:222HHLLT4:2:1154:21231:18398 113     HOFI_CYP2J19    116     60      151M    =       116     -151GACCCAGGAATTTCCCTCCAGGGCCGCAGCTCTTTCCTCTCGTGGGAACCTTTGTGGACTTTAAGCAGCCCCTCCATCTTGCACTGCAGAAGCTTACGGGTCGGTACGGGAACATCTTCAGCGTGCAGTTTGGAAGTCTGACATTCGTGGT     9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     NH:i:1  AS:i:146    NM:i:1
LH00260:42:222HHLLT4:2:1154:21231:18398 177     HOFI_CYP2J19    116     60      151M    =       116     -151GACCCAGGAATTTCCCTCCAGGGCCGCAGCTCTTTCCTCTCGTGGGAACCTTTGTGGACTTTAAGCAGCCCCTCCATCTTGCACTGCAGAAGCTTACGGGTCGGTACGGGAACATCTTCAGCGTGCAGTTTGGAAGTCTGACATTCGTGGT     9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     NH:i:1  AS:i:146    NM:i:1
```
If we want to visualize these reads we can convert the sam to a bam file and load it into IGV. 

```
samtools view -b -S -T CYP2J19.fasta HEK_CYP2J19_blast.sam -o HEK_CYP2J19_blast.bam

samtools sort HEK_CYP2J19_blast.bam -o HEK_CYP2J19_blast_sort.bam

samtools index HEK_CYP2J19_blast_sort.bam
```

* [blast_bam.sh](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/blast_bam.sh)
* [blast_bam.sbatch](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/blast_bam.sbatch)

Now download `CYP2J19.fasta`, `CYP2J19.fasta.fai`, `HEK_CYP2J19_blast_sort.bam`, and `HEK_CYP2J19_blast_sort.bam.bai` to you PC and load them into IGV using `CYP2J19.fasta` as the "genome". 

![](https://github.com/mbtoomey/genome_biology_FA24/blob/main/Lessons/scripts/annotation_1.png)

Now we can see where the reads are mapping within this particular transcript and pick out SNPs. You might notice that is redundant with the mapping approaches we have discussed elsewhere. However, magicblast may perform better with error prone sequencing (i.e. nanopore reads) and may offer imporved intron detection [(Boratyn et al. 2019)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2996-x).





