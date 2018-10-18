# FACETS
This repository contains an implementation of [FACETS](https://github.com/mskcc/facets), an algorithm to implement fraction and copy number estimate from tumor/normal sequencing, for FireCloud. FACETS is run on tumor-normal pairs to determine the purity and ploidy of a tumor sample, as well as infer allele-specific copy number. 

The first step of the algorithm performs a pileup to obtain coverage information at many sites of common single nucleotide variants in both the normal bam and tumor bam. While resources such as [dbSNP](https://www.ncbi.nlm.nih.gov/projects/SNP/) would work, the [generate hetsites folder](https://github.com/vanallenlab/facets/tree/master/generate_hetsites) features an ipython notebook to generate a bed file of common sites from [ExAC](http://exac.broadinstitute.org/) and instruction on how to generate a VCF of just those sites. This output, exac_hets.vcf, is [publicly available on Google Cloud](https://console.cloud.google.com/storage/browser/fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/facets).

**Please read the manuscript and the [userguide](https://github.com/mskcc/facets/blob/master/vignettes/FACETS.pdf) before using**. 

Github: [mskcc/facets](https://github.com/mskcc/facets)  
DOI: [https://doi.org/10.1093/nar/gkw520](https://academic.oup.com/nar/article/44/16/e131/2460163)  
Citation: [Ronglai Shen, Venkatraman E. Seshan; FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing, Nucleic Acids Research, Volume 44, Issue 16, 19 September 2016, Pages e131, https://doi.org/10.1093/nar/gkw520](https://academic.oup.com/nar/article/44/16/e131/2460163)

Docker image: [vanallenlab/facets](https://hub.docker.com/r/vanallenlab/facets/)  
FireCloud method: [breardon/facets](https://portal.firecloud.org/#methods/breardon/facets/)
