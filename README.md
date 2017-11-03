![alt text](http://empede.co.uk/imgrepos/OMGene_head.png? "OMGene logo")

What does OMGene do?
==========
OMGene (Optimise My Gene) aims to improve simultaneously sets of gene models on an orthogroup-by-orthogroup basis. It leverages genetic information from multiple orthologs across multiple genomes to find the most evolutionary consistent set of gene models for a given gene.

Given an input set of gene models and a set of genomes, it will seek to find the gene model variants for each gene which collectively optimise the coherence of the multiple sequence alignment, thus maximising the likelihood that a gene model is correct in the absence of extrinsic data. OMGene requires only a set of single-entry GTF gene model files and a set of input FASTA genome files.

**Github link:** https://github.com/mdunne/OMGene

![alt text](http://empede.co.uk/imgrepos/OMGene_examples.png "OMGene examples")


Usage
=====

### Default usage

OMGene runs as a single command that takes as input a tabulated text file containing locations of genomes and GTF files in the following format:

```
#gtf                                  genome
/path/to/gtf_species1_gene1.gtf       /path/to/genome_species1.fasta
/path/to/gtf_species2_gene1.gtf       /path/to/genome_species2.fasta
/path/to/gtf_species3_gene1.gtf       /path/to/genome_species3.fasta
etc.
```

Each GTF file should contain the coordinates of a *single gene*. OMGene is then run using:

`python omgene.py -i path/to/gene_info_file.tdv -o output_folder`

If no output folder is specified then OMGene will create one with a generic name. Currently OMGene is only configured to run on a single core.

### Example files

An example file is given in the "samples" directory included in this repository. This contains actual published gene models and truncated FASTA genome sections. To run the example, navigate to the samples directory and execute the following command:

`python ../omgene.py -i in.tdv -o out`

Once this has run, the "out" directory will contain the results. To compare the original and new results, compare the files working/or.aln and results/all.aln. To view this quickly in the terminal, you may wish to use the [Alan alignment viewer](https://github.com/mpdunne/alan). If you have alan installed, you can type

`cat working/or.aln results/all.aln > compare.aln; alan compare.aln`

### Alternate start and splice sites

You may wish to use non-canonical nucleotide combinations for start and splice sites, for example GU--GG splicing or CUG start codons. These can be specified using the `-do`, `-ac`, and `-sc` options for donor sites, acceptor sites, and start codons respectively. Site lists must be comma-separated.If you wish to include the canonical codon or splice site in the search, you must also specify it here. For example:

`python ../omgene.py -i in.tdv -o out -do "gc,gu" -ac "gg,ag" -sc "aug,cug"`.

Cysteine and uracil (C and U) will be interpreted identically. If these options are omitted, the default will be used: GU,GC donor sites, AG acceptor sites, and ATG start codons. The more options that are inputted here, the more variants will be probed by OMGene, potentially considerably increasing its runtime. Therefore we recommend some degree of restraint here when choosing your start and splice sites.


### Search margins

By default, OMGene searches the region 600bp either side of the termini of each inputted gene. This buffer margin can be expanded or contracted as desired using the -b option, e.g.

`python omgene.py -i path/to/gene_info_file.tdv -o output_folder -b 1000`

### Serving suggestion

It is recommended that orthogroups be obtained using [OrthoFinder](https://github.com/davidemms/OrthoFinder), implemented on primary transcripts only.

The results from OrthoFinder are outputted as a large, tab-delimited table of gene IDs. You will need to use these IDs to extract genomic GTF coordinates from your species GTF file, and produce a TDV input file as described above.

### Reference species

If you do not wish to change all of the gene models in your species group, you may divide them into *target* and *reference* species. Target species are specified in a file given by the `-i` option, and reference species are specified in the file given by the `-r` option.

`python omgene.py -i path/to/target_species.tdv -r path/to/reference_species.tdv -o output_folder`

The `target_species.tdv` file should *only* contain the species whose genomes you would like to search, and the `reference_species.tdv` file should contain *only* species which you *do not* wish to search.

### OMGene likes clean FASTA files.

The current iteration of OMGene requires that chromosome names in the genome FASTA files consist of unspaced IDs only (i.e. no description lines or information other than the name). This will be fixed in future versions. In the meantime, you may wish simply to modify the names of your chromosome entries using:

```
sed -r "s/>([^ ]+*) .*/>\1/g" genome.fa > genome_clean.fa
```

Obtaining GTFs from GFF3 files
==============================
OMGene uses GTF files for both its input and output, due to the superior uniformity of GTF naming and attribute conventions. To convert files from GFF3 to GTF format, we recommend using the simple tool fml_gff3togtf, from the Galaxy Tool Shed. This can be found at https://toolshed.g2.bx.psu.edu/repository?repository_id=afcb6456d8e300ed, and is implemented using python:

```
python gff_to_gtf.py infile.gff3 > outfile.gtf
```

Users should note that this tool truncates chromsome names to 15 characters. If this is going to be an issue, a wrapper for this script can be found in the utils directory in this repository (https://github.com/mpdunne/orthofiller/blob/master/utils/gff_to_gtf_safe.py). The above Galaxy tool should be downloaded first, and the path to its directory should be included in the appropriate place at the top of the `gff_to_gtf_safe.py` file. The full script can then be run as

```
python gff_to_gtf_safe.py infile.gff3 outfile.gtf
```

Note that, in order to function properly, the above conversion script requires that entries in the GFF3 input file are well-formed: that is they contina gene, mRNA, CDS, and exon entries for each gene. Ideally ensure that each GFF3 entry has each of these attributes before proceeding. Alternatively, if you simply wish to remove incomplete entries from your GFF3 file, you can use the `clean_gff.py` script, also included in the utils directory of this repository. The usage for this script is:

```
python clean_gff.py infile.gff3 infile_clean.gff3
```

Output File Format
==================
OMGene output can be found in the `results` folder of the specified output directory. For each inputted gene, each of a CDS and amino acid FASTA file are produced, as well as a GTF file containing the genomic coordinates of the relevant gene.

Installing Dependencies
=======================
OMGene is written to run on linux and requires the following to be installed and in the system path:

1. Python 2.7 together with the scipy, numpy, pyfaidx, networkx and BioPython libraries 

2. BedTools Suite

3. Exonerate

### python and scipy

Up-to-date and clear instructions are provided here: http://www.scipy.org/install.html, be sure to chose a version using python 2.7. As websites can change, an alternative is to search online for "install scipy".

### BedTools

The BedTools suite can be downloaded from [https://github.com/arq5x/bedtools2/releases](https://github.com/arq5x/bedtools2/releases)

### Exonerate

Exonerate can be downloaded from the EBI website, https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate


Other useful software
=====================

### OrthoFiller

OrthoFiller is a pipeline designed to help "fill in" gaps in the proteomes of a set of species, by searching for absent members of orthogroups.

https://github.com/mpdunne/orthofiller

### Alan

Alan is an in-terminal command-line tool for viewing alignments without the need for a GUI.

https://github.com/mpdunne/alan

