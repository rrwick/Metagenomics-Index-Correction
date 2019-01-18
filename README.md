# GTDB taxonomy scripts

This repo contains scripts used to prepare and analyse Centrifuge classifications using the [GTDB database](http://gtdb.ecogenomic.org/). They were used in the preparation of the __TITLE-OF-PAPER__ manuscript.



## tax_from_gtdb.py

This script generates an [NCBI-style](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt) taxonomy using the GTDB database and other files necessary for building a [Centrifuge](https://ccb.jhu.edu/software/centrifuge) or Kraken database using GTDB.



## count_classifications.py

This script takes a Centrifuge SAM file as input and produces two outputs:
  * (to stdout) a single-row table summarising the read set
  * (to stderr) a one-row-per-read table describing each read

For reads with multiple classifications to different tax IDs, this script will find the LCA of the tax IDs to give a final classification. It also categories each read based on the taxonomic rank of its final classification (unclassified, root, domain, phylum, class, order, family, genus or species).

The NCBI taxonomy contains many additional and in-between ranks (e.g. subspecies, superfamily). If the read's final classification falls into one of these then the script will give it the rank of the first ancestor with a standard rank. E.g. subspecies -> species, superfamily -> order.



## dereplicate_assemblies.py

This script uses Mash-based clustering to dereplicate GTDB assemblies. This allows for a major reduction in the number of assemblies (making for faster and smaller databases) which still retaining much of the genomic diversity.



## read_set_n_count.py

This is a simple script written to deal with the fact that some [HMP](https://hmpdacc.org/) read sets are full of `N`-heavy reads. It takes one or more read files as input and outputs a table which shows information on the number of reads which are mostly `N`.



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
