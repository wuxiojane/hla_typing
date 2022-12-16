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

The input file format for this process can be fastq and bam. The bam file that mapped to genome is not recommended. Due to the complexity of HLA gene regions, many reads cannot be successfully compared to HLA gene regions. Recommanded the bam files which mapped with the IPD-IMGT/HLA database.
