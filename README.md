# binf6110-assignment-3

## Introduction
One study conducted in 2019 by De Filippis et al., which investigated human gut microbomes of Italians using shotgun metagenomics, found that differing dietary habits can actually affect the presence of certain microbiota, especially for one prevalent species, *Prevotella copri*.  They found that this species revealed different strains at different locations based on the types of foods consumed often, such as fiber-rich or omnivorous, revealing how diet can affect diversity at even the finest levels of taxonomy (De Filippis et al., 2019).  Interestingly enough, it also appeared to suggest how diets between different types of global populations, and in turn, the composition of gut microbiota, can further affect the ability of biological processes to consume certain types of food, most notably in non-Western countries like Italy (De Filippis et al., 2019).

Studies such as these can most useful for explaining the essence of these connections between microbial taxonomy and lifestyle in human health/disease as a whole and across populations (De Filippis et al., 2019).  However, it is also important that the composition and abundance of other taxa within the diversity should be investigated more thoroughly, to understand exactly how they collectively support the necessary biological functions for food consumption and metabolism.  Hence, this taxonomic analysis seeks to expand on the significance of these connections, by applying different classification and diversity pipelines to compare how the abundances and diversity of certain microbiotic taxa in the human gut metagenomes differ across individuals with varying dietary patterns.  It will also attempt to explain how the presence and/or dominance of certain species across diets support digestive and metabolic functions and processes most associated with them.


## Methods
Raw short-read data for human gut microbiomes was randomly selected from samples used in a shotgun metagenomics study, with Illumina sequencing, from 74 Italians based on dietary habits, found on the National Center for Biotechnology Information (NCBI) database (University of Naples Federico II, 2017).  A total of six sets of paired-end samples were obtained, with three sets belonging to individuals with vegan diets and the other three belonging to individuals with omnivorous diets (University of Naples Federico II, 2017).  Quality control was conducted on each individual paired-end sample using FastQC (v0.12.1) on 12 threads per sample, to check for sources of contamination, which were then compiled into overall reports with MultiQC (v1.33) (Andrews, 2010, Ewels et al., 2016).  Adapter content was then trimmed off each paired-end sample using *fastp*, before repeating quality control reports a second time and compressing them (.gz files) (Chen et al., 2018).  Taxonomic classification for each set of paired-end samples (.gz) was done using k-mer classification in Kraken2 (v2.1.6), using a confidence level of 0.15 to reduce the likelihood of false positives (Wood et al., 2019).  These computations were run remotely in the Narval Compute Canada cluster at Digital Research Alliance of Canada (DRAC), with the most recent Kraken2 Standard database as of February 26, 2026 (Langmead, 2026).  Finally, taxonomic abundances for each sample were re-estimated at species level from classification genus reads, using Bayesian Reestimation of Abundance with KrakEN (Bracken v3.0.1) (Lu et al., 2017).  These Bracken species reports were subsequently converted into BIOM files using *kraken-biom* (v1.2.0) for further taxonomic analysis (Dabdoub, 2016).

Continuing in R (v4.5.2), all BIOM reports of abundance tables from Bracken were imported into the working environment with the *phyloseq* package (v1.54.2) and the *biomformat* package (v1.38.3) (McMurdie & Holmes, 2013, McMurdie & Paulson, 2026).  Relative abundance of taxa between dietary groups and samples were compared and visualized with barplots using *phyloseq* (v1.54.2) to explore the diversity of classifications and determine the most prominent microbial species in the metagenomes (McMurdie & Holmes, 2013).  Alpha and beta diversity measures were both calculated using functions from the *vegan* package (v2.7-3) (Oksanen et al., 2026).  For alpha diversity, the Shannon index was calculated to understand and visualize the ovarall taxa richness, while the Simpson index was calculated to further evaluate taxa evenness and/or dominance between groups and/or samples.  Meanwhile, beta diversity was compared between vegan and omnivore diet groups using Jaccard and Bray-Curtis distance measures for similarities and dissimilarities in microbiomes.  Lastly, to determine the key differences in taxa, ANCOM-BC was used to conduct differential abundance analysis between the two groups (). 

## Results

## Discussion

## References
Andrews, S. (2010). FastQC. A quality control tool for high throughput sequence data. *Babraham Bioinformatics*. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, *34*(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

Dabdoub, S. M. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2). https://github.com/smdabdoub/kraken-biom

De Filippis, F., Pasolli, E., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, *25*(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, *32*(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

Langmead, B. (2026). Index zone - Kraken 2 / Bracken Refseq indexes. https://benlangmead.github.io/aws-indexes/k2

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, *3*, e104. https://doi.org/10.7717/peerj-cs.104

McMurdie, P. J. & Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLOS One*, *8*(4), e61217. https://doi.org/10.1371/journal.pone.0061217

McMurdie P, Paulson J (2026). biomformat: An interface package for the BIOM file format [R package version 1.38.3]. https://github.com/joey711/biomformat/
  
Oksanen, J., Simpson, G., Blanchet, F., Kindt, R., Legendre, P., Minchin, P., O'Hara, R., Solymos, P., Stevens, M., Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B., Borcard, D., Borman, T., Carvalho, G., Chirico, M., De Caceres, M., Durand, S., Evangelista, H., FitzJohn, R., Friendly, M., Furneaux, B., Hannigan, G., Hill, M., Lahti, L., Martino, C., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., Smith, T., Stier, A., Ter Braak, C., & Weedon, J. (2026). vegan: Community Ecology Package [version 2.7-3]. https://doi.org/10.32614/CRAN.package.vegan

University of Naples Federico II. (2017). SRX4967482: gut metagenome of vegan subject from Bari (Italy), subject ID 01BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967482[accn]

University of Naples Federico II. (2017). SRX4967479: gut metagenome of vegan subject from Parma (Italy), subject ID 26PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967479[accn]

University of Naples Federico II. (2017). SRX4967478: gut metagenome of omnivore subject from Bari (Italy), subject ID 03BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967478[accn]

University of Naples Federico II. (2017). SRX4967471: gut metagenome of vegan subject from Turin (Italy), subject ID VOV56. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967471[accn]

University of Naples Federico II. (2017). SRX4967459: gut metagenome of omnivore subject from Turin (Italy), subject ID VOV77. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967459[accn]

University of Naples Federico II. (2017). SRX3463060: gut metagenome of omnivore subject from Parma (Italy), subject ID 37PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX3463060[accn]

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, *20*, 257. https://doi.org/10.1186/s13059-019-1891-0
