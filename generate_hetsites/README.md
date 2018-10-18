# Generate hetsites for pileup
We use SNPs from ExAC that are more common than allele frequencies of 10^-2 for chromosomes 1-X and more common than 10^-5 for Y. 

Create exac_hets.vcf with GATK 3.8
`java -jar ~/software/gatk/GenomeAnalysisTK-3.8.jar -T SelectVariants -R ~/storage/hg19/Homo_sapiens_assembly19.fasta -V ~/storage/exac/ExAC.r1.sites.vep.vcf.gz -o output.vcf -L exac_hets.bed`
