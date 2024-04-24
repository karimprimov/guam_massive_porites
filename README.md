In this script, code is included for demultiplexing, trimming, generating multi-locus genotypes, and running population genomic analyses. This code includes examples of each step in the Linux bioinformatic and R pipeline. This code was produced to study the population genomics of a morphologically convergent coral species group called massive Porites. Found throughout the Indo-Pacific, they occur in different reef environments, and in this study, we find that different species,not populations of the same species, predominantly occupy river deltas and adjacent (< 100M) fore reefs on Guam respectively, and bleaching less prevalent in river deltas. 

In this analysis, we use the following programs to curate our datasets, and subsequently analyse our data in both Linux and RStudio:

1. Custom demultiplexing script from Hannah Weigand
2. Bowtie2 (building reference database & aligning our sample files to the reference genome)
3. Samtools to convert sam to bam files, and then sort & index them
4. Generate multi-locus genotypes calls using STACKS (genotypes calls) & ANGSD (genotypes likelihoods)
5. Phylogenetic analyses using the CIPRES online phylogenetic tree website (input data was a STACKS-generated phylip file)
7. Population genomic analyses using both STACKS & ANGSD
