# hla_typing
After comparing the following softwares: 
* OptiType
* Phlat
* HLA-HD
* HISAT-genotype
* hlavbseq v2
* hla-la  

It is concluded that the **HLA-HD** has a better performance.

## Dependencies 
The HLA typing process of WGS data was built based on HLA-HD software, which depended on **bowtie2**, **samtools** and **HLA-HD**.

The input file format for this process can be fastq and bam. The bam file for genomic comparison is not recommended. Due to the complexity of HLA gene regions, many reads cannot be successfully compared to HLA gene regions. bam files can be used against the IPD-IMGT/HLA database.
