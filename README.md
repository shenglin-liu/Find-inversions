# Find-inversions
Find inversions in the genome using high-density SNP data of multiple populations

# Description
It is an R code developed to detect inversions in the genome of a diploid species. It needs high density SNP data (RAD or whole genome sequences) of multiple populations. The inversions need to have different frequencies across the populations in order to be detected. For this method to apply, a reference genome is required, or at least the relative positions of the SNPs should be available.

It takes a VCF file as input, and it also needs a popmap file to provide population affinity of the individuals in the VCF file. The popmap file should contain two columns delimited by tabs. The first column enlists the individuals in the VCF file in the same order. The second column gives the population name for each individual.

# Principle
In an inversion, the recombination between the inverted and the non-inverted haplotypes is suppressed. As a result, many SNPs in this region should share the same allele frequency in a population, especially if the two haplotypes have diverged for a long time. Consequently, in a Manhattan plot of FST between two populations exhibiting different freqeuncies of the two haplotypes, we should see a rectangular peak of FST for neighboring SNPs in the inverted region.

# Algorithm
The algorithm is based on FST calculation. But the roles of populations and SNPs in the calculation are reversed. Conventionally, FST is calculated between two populations across many SNPs. But in this algorithm FST is calculated between SNPs across many populations. Thus the SNPs sharing the same allele frequency, i.e., the rectangular peaks, will have consistently low FST across populations.

# Performance
The accuracy of the analysis is mainly influenced by the following factors.

Number of populations: The more the better; minimum is three.

Heterogeneity of the frequency of the inverted haplotype across the populations: The higher the better.

Accuracy of the allele frequency estimates: The higher the better; it requires decent sample sizes and high quality of SNP calling.

# Instructions
Download the R script to the working directory of your R session. Modify the "Parameters" section according to your need. Run the script using source() function. Here are some guidelines about how to set the parameters.

f.vcf: Name of the VCF file.

f.popmap: Name of the popmap file.

referencePop: Name of the population used for defining the non-inverted haplotype; it is arbitrary.

n.digit: Number of digits to be kept in calculating allele frequencies; VERY IMPORTANT, but the default (4) is fine unless you have very big sample sizes.

n.nonZeroPop: A SNP must be polymorphic in at least this number of populations, otherwise it will be filtered out; SNPs showing polymorphism in very few populations can easily cause false positive.

lowFreq: VERY IMPORTANT!!! SNPs with overall frequency lower than this will be removed.

l.window: Sliding window size; default is fine unless the genome assembly is very fragmented.

l.step: Sliding step; default is fine unless the genome assembly is very fragmented.

fst.threshold: VERY IMPORTANT!!! FST threshold; FST higher than this will be removed.

minRange: Minimum size of an inversion.

minInter: Minimum distance between neighboring informative SNPs (those giving rectangular peaks); if two SNPs have less than this distance, only one will be used for inerence, but both will be reported in the result.

minSnps: Minimum number of informative SNPs for an inversion; should be no less than 3.

f.output: The name of the output file.

# Output
The output file enlists the detected inversions in blocks separated by two empty lines. Each block contains the inversions detected in a chromosome. Each inversion is specified by starting with the starting position and ending positing of the inversion. A table follows listing the informative SNPs (those giving rectangular peaks) in the inversion and giving the frequency of the inverted copy.

# Citation
Liu S, Ferchaud AL, Groenkjaer P, Nygaard R, Hansen MM (2017) Genomic parallelism and lack thereof: Determinants inferred from contrasting systems of three-spine sticklebacks. Submitted.
