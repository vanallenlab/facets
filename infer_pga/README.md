# Infer percent genome altered from FACETS

Infers the percent genome altered using FACETS total copy number. Fundamentally, the script takes all the segments that do not match diploid copy number (2 for autosomes and 1 for sex chromosomes) and computes their sizes, then divides by the total size of all the segments.

The output is written to `facets_fraction_genome_CNA_altered.txt`





