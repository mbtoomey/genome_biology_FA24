# BLAST

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

To create a database of RNA sequences to search, I download all of the RNAs from the NCBI human refernce genome - [GCF_000001405.40_GRCh38.p14_rna.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz)

