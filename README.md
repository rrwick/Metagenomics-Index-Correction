# Metagenomics Index Correction

This repository contains scripts used to prepare, compare and analyse metagenomic classifications using custom index databases, either based on default [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) or [GTDB](http://gtdb.ecogenomic.org/) taxonomic systems. These scripts were used in the preparation of the Méric, Wick et al. (2019) manuscript entitled: [_Correcting index databases improves metagenomic studies_](https://www.biorxiv.org/content/10.1101/712166v1).

Our custom index databases are hosted on figshare:<br>
[monash.figshare.com/projects/Metagenomics_Index_Correction/65534](https://monash.figshare.com/projects/Metagenomics_Index_Correction/65534)


## Table of contents

* [Using custom indices](#using-custom-indices)
  * [Centrifuge](#centrifuge)
  * [Kraken2](#kraken2)
* [Scripts for analysis/index creation used in the manuscript](#scripts-for-analysisindex-creation-used-in-the-manuscript)
  * [Requirements](#requirements)
  * [Custom classification indices using GTDB definitions](#custom-classification-indices-using-gtdb-definitions)
  * [Taxonomic rank counts from Centrifuge output](#taxonomic-rank-counts-from-centrifuge-output)
  * [Custom dereplication thresholds for GTDB assemblies](#custom-dereplication-thresholds-for-gtdb-assemblies)
  * [Counting N-heavy reads in FASTQ sequencing reads](#counting-n-heavy-reads-in-fastq-sequencing-reads)
  * [Finding reads unclassified in one index but not another](#finding-reads-unclassified-in-one-index-but-not-another)
* [License](#license)

## Using custom indices

### Centrifuge

We used [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) in our manuscript because it uses less RAM than alternative tools. Centrifuge indices consist of a few files named `*.cf`.

Please consult the [Centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml) for full instructions on its usage, but in brief the steps are:

1. Download a pre-compiled index (from our [our figshare project](https://monash.figshare.com/projects/Metagenomics_Index_Correction/65534)).
    * Latest (recommended) index: [`GTDB_r89_54k`](https://monash.figshare.com/articles/GTDB_r89_54k/8956970) 
    * Older indices described in the manuscript: [`GTDB_r86_46k & NCBI_r88_Human17k`](https://monash.figshare.com/projects/Metagenomics_Index_Correction/65534) 
2. Run Centrifuge on a paired-end Illumina read set:<br>
`centrifuge -p 16 -x index_dir -1 *_1.fq.gz -2 *_2.fq.gz`<br>
`-S classifications --report-file report.tsv`
3. Generate a Kraken-style summary report:<br>
`centrifuge-kreport -x index_dir classifications > kreport.txt`


### Kraken2

[Kraken2](https://ccb.jhu.edu/software/kraken2/) is an alternative metagenomics classification tool. It is faster than Centrifuge but may use more RAM. Its indices consist of a few `*.k2d` files. Please consult the [Kraken2 manual](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) for usage instructions. As part of this work, we generated Kraken2 (and Kraken1) custom indices too, which can downloaded from [our figshare project](https://monash.figshare.com/projects/Metagenomics_Index_Correction/65534).




## Scripts for analysis/index creation used in the manuscript

### Requirements

To run these scripts you will need:
* Python 3.5 (or later)
* [Biopython](https://biopython.org/)
* [Mash](https://github.com/marbl/Mash) installed (able to be run from the command line)


### Custom classification indices using GTDB definitions

The `tax_from_gtdb.py` script generates [NCBI-style](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt) taxonomy files using the [GTDB](http://gtdb.ecogenomic.org) definitions. It can also prepare the sequence files necessary for building a Centrifuge and/or Kraken2 database.

Usage – just make the taxonomy files:
```
tax_from_gtdb.py --gtdb taxonomy.tsv --nodes nodes.dmp --names names.dmp
```

Usage – prepare for [building a Centrifuge database](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database):
```
tax_from_gtdb.py --gtdb taxonomy.tsv --assemblies genomes --nodes ex.tree --names ex.name --conversion ex.conv --cat_fasta ex.fa
```

Usage – prepare for [building a Kraken2 database](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases):
```
tax_from_gtdb.py --gtdb taxonomy.tsv --assemblies genomes --nodes nodes.dmp --names names.dmp --kraken_dir kraken_genomes
```

Input:
* `--gtdb`: taxonomy file obtained from the [downloads section of the GTDB website](http://gtdb.ecogenomic.org/downloads). Required. GTDB provides separate taxonomy files for bacteria and archaea, but these can be concatenated together to make a single taxonomy file for this script.
* `--genomes`: a directory containing all the reference genomes sequences in FASTA format (can be gzipped). Required if you are preparing for a Centrifuge/Kraken2 database build.

Output:
* `--nodes`: filepath to save the nodes file (equivalent to NCBI's `nodes.dmp`).
* `--names`: filepath to save the names file (equivalent to NCBI's `names.dmp`).
* `--conversion`: filepath to save the [conversion table for Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database).
* `--cat_fasta`: filepath to save the [concatenated reference seq file for Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database)
* `--kraken_dir`: a directory to save the [Kraken2-ready reference FASTAs](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases).




### Taxonomic rank counts from Centrifuge output

The `count_classifications.py` script allows the breakdown of read classification to different taxonomic ranks from a raw classification output from Centrifuge. 

For reads with multiple classifications to different tax IDs, this script will find the LCA of the tax IDs to give a final classification. It also categories each read based on the taxonomic rank of its final classification (unclassified, root, domain, phylum, class, order, family, genus or species).

The NCBI taxonomy contains many additional and in-between ranks (e.g. subspecies, superfamily). If the read's final classification falls into one of these then the script will give it the rank of the first ancestor with a standard rank. E.g. subspecies -> species, superfamily -> order.

Usage:
```
count_classifications.py --centrifuge classifications --tree nodes.file --prefix out_prefix 1> read_set_summary 2> per_read_summary
```

Input:
* `--centrifuge`: the output of Centrifuge in tab-delimited format
* `--tree`: the taxonomy file used to build the Centrifuge index

Output:
  * `--prefix`: prefix of the tsv files this script will make
  * `read_set_summary`: a single-row table summarising the read set (to stdout)
  * `per_read_summary`: a one-row-per-read table describing each read (to stderr)


### Custom dereplication thresholds for GTDB assemblies

"Dereplication" refers to the threshold-based selection of representative reference genomes for phylogenetically-similar clusters, and has been used in the GTDB study to provide a reduced, less-redundant list of 28941 bacterial and archaeal reference genomes that are representative of similarity clusters on a phylogenetic tree. By default, the "dereplicated" list proposed on the [downloads section of the GTDB website](http://gtdb.ecogenomic.org/downloads) contains (release 86; September 2018) 27372 bacterial and 1569 archaeal genomes. As explained in the [GTDB publication](https://www.nature.com/articles/nbt.4229), the dereplication was made according to various criteria, but mostly defining two genomes as being "replicates" when their Mash distance was ≤0.05 (∼ANI of 95%).

Here, our `dereplicate_assemblies.py` script uses Mash-based clustering to dereplicate GTDB assemblies using a user-defined Mash threshold. This allows for a reduction in the number of assemblies (making for faster and smaller Centrifuge index databases), while still retaining a variable (user-defined) amount of genomic diversity.

Usage:
```
dereplicate_assemblies.py --threads 32 --threshold 0.005 all_assemblies_dir derep_assemblies_dir bac_and_arc_taxonomy_r86.tsv
```

Parameters:
* `--threads`: number of CPU threads
* `--threshold`: Mash distance threshold. The GTDB_r86_46k index from the associated publication used 0.005.

Input:
* `all_assemblies_dir`: path to all assemblies are located, typically all assemblies from GenBank/RefSeq.
* `bac_and_arc_taxonomy_r86.tsv`: GTDB taxonomy for all bacterial genomes from the [GTDB website](http://gtdb.ecogenomic.org/downloads)

Output:
* `derep_assemblies_dir`: output path


### Counting N-heavy reads in FASTQ sequencing reads

The simple `read_set_n_count.py` script has been written to deal with the fact that some [HMP](https://hmpdacc.org/) read sets are full of `N`-heavy reads. It takes one or more read files as input and outputs a table which shows information on the number of reads which are mostly `N`.

Usage:
```
read_set_n_count.py file.fastq > counts
```
Input:
* `file.fastq`: read file to be counted

Output: 
* `counts`: a tab-delimited file with the following four columns: filename, number of reads, number of N reads, percentage of N reads.<br>
e.g. `file.fastq    44880075    44382676    98.892%`


### Finding reads unclassified in one index but not another

The `find_unclassified.py` script can compare Centrifuge classifications between two indices and identify reads that are unclassified using one index but not another. In order to make sense, this script must be run using the outcomes of two classifications of the same sample using two different indices. This script was used in our analysis to reclassify using the `nt` database reads that had remained unclassified using the GTDB_r86_46k index.

Usage:
```
find_unclassified.py index1_centrifuge_output index2_centrifuge_output | grep -v unclassified > output
```
Input:
* `index1_centrifuge_output` and `index2_centrifuge_output`: the two classification files to be compared

Output:
* `output`: classification in the same format as the regular Centrifuge output (to stdout)



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
