# Metagenomics Index Correction
This repository contains scripts used to prepare and analyse Centrifuge classifications using custom index databases, either based on default NCBI or [GTDB](http://gtdb.ecogenomic.org/) taxonomic systems. These scripts were used in the preparation of the Méric, Wick et al. (2019) manuscript entitled: "XXXXXX". Brief instructions are provided below each scripts.



## Creating a Centrifuge custom index using GTDB taxonomic definitions
The `tax_from_gtdb.py` script generates an [NCBI-style](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt) taxonomy using the GTDB definitions (from the GTDB website) and other input files necessary for building a [Centrifuge](https://ccb.jhu.edu/software/centrifuge) or Kraken index database using GTDB.

Usage:

`./tax_from_gtdb.py --gtdb taxonomy.tsv --assemblies PATH --nodes nodes.file --names names.file --conversion conv.file --cat_fasta fasta.file`

With: 

`taxonomy.tsv`: taxonomy file obtained from the [downloads section of the GTDB website](http://gtdb.ecogenomic.org/downloads).

`PATH`: the pathname to all reference genomes sequences in compressed or uncompressed fasta format.

`nodes.file`: name and path of the output nodes file required to build an index with Centrifuge

`names.file`: name and path of the output names file required to build an index with Centrifuge

`conv.file`: name and path of the output conversion file required to build an index with Centrifuge

`fasta.file`: name and path of the output fasta file required to build an index with Centrifuge



## Taxonomic rank count from Centrifuge classification output
The `count_classifications.py` script allows the breakdown of read classification to different taxonomic ranks from a raw classification output from Centrifuge. It takes a Centrifuge SAM file as input and produces two outputs:
  * (to stdout) a single-row table summarising the read set
  * (to stderr) a one-row-per-read table describing each read

For reads with multiple classifications to different tax IDs, this script will find the LCA of the tax IDs to give a final classification. It also categories each read based on the taxonomic rank of its final classification (unclassified, root, domain, phylum, class, order, family, genus or species).

The NCBI taxonomy contains many additional and in-between ranks (e.g. subspecies, superfamily). If the read's final classification falls into one of these then the script will give it the rank of the first ancestor with a standard rank. E.g. subspecies -> species, superfamily -> order.

Usage:

`./count_classifications.py --sam file.sam --tree nodes.file --prefix PREFIX`

With:

`file.sam`: the output of Centrifuge in SAM format

`nodes.file`: the file used to build the Centrifuge index containing taxonomic information

`PREFIX`: prefix of the output files of this script



## Creating custom "dereplication" thresholds for assemblies used to build a GTDB-based Centrifuge index

"Dereplication" refers to the threshold-based selection of representative reference genomes for phylogenetically-similar clusters, and has been used in the GTDB study to provide a reduced, less-redundant list of 28941 bacterial and archaeal reference genomes that are representative of similarity clusters on a phylogenetic tree. By default, the "dereplicated" list proposed on the Downloads section of the [GTDB website](http://gtdb.ecogenomic.org/downloads) contains (release 86; September 2018) 27372 bacterial and 1569 archaeal genomes. As explained in the [GTDB publication](https://www.nature.com/articles/nbt.4229), the dereplication was made according to various criteria, but mostly defining two genomes as being "replicates" when their Mash distance was ≤0.05 (∼ANI of 95%).

Here, our `dereplicate_assemblies.py` script uses Mash-based clustering to dereplicate GTDB assemblies using a user-defined Mash threshold. This allows for a major reduction in the number of assemblies (making for faster and smaller Centrifuge index databases), while still retaining a variable (user-defined) amount of genomic diversity.

Usage:

`./dereplicate_assemblies.py --threads N --threshold 0.005 all_assemblies_dir derep_assemblies_dir bac_and_arc_taxonomy_r86.tsv`

With:

`--threads N`: number of threads, e.g. 32.

`--threshold`: Mash distance threshold defining. The GTDB_r86_46k index from the associated publication used 0.005.

`all_assemblies_dir`: path to all assemblies are located, typically all assemblies from GenBank/RefSeq.

`derep_assemblies_dir`: output path

`bac_and_arc_taxonomy_r86.tsv`: GTDB taxonomy for all bacterial genomes from the [GTDB website](http://gtdb.ecogenomic.org/downloads)

## Counting "N-heavy" reads in FASTQ sequencing reads.

The simple 'read_set_n_count.py' script has been written to deal with the fact that some [HMP](https://hmpdacc.org/) read sets are full of `N`-heavy reads. It takes one or more read files as input and outputs a table which shows information on the number of reads which are mostly `N`.

Usage:

`./read_set_n_count.py file.fastq`

Output: 

four tab-delimited columns:
          1) filename
          2) number of reads
          3) number of N reads
          4) percentage of N reads

e.g. `file.fastq    44880075    44382676    98.892%`


## Finding reads unclassified in index1 but not in index2 (for the same sample)

To compare Centrifuge classifications between two indices, and identify reads in one particular sample that are unclassified (or another outcome) in one index compared to another, the `find_unclassified.py` script can be used. In order to make sense, this script must be run using the outcomes of two classifications of the same sample using two different indices. This script was used in our analysis to reclassify using the `nt` database reads that had remained unclassified using the GTDB_r86_46k index.

Usage:

`./find_unclassified.py sample_index1.sam sample_index2.sam | grep -v "unclassified" > output.sam`

Output:

`output.sam`: similar SAM format as the regular Centrifuge output.


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
