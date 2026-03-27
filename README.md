# binf6110-assignment-3

## Introduction
One study conducted in 2019 by De Filippis et al., which investigated human gut microbomes of Italians using shotgun metagenomics, found that differing dietary habits can actually affect the presence of certain microbiota, especially for one prevalent species, *Prevotella copri*.  They found that this species revealed different strains at different locations based on the types of foods consumed often, such as fiber-rich or omnivorous, revealing how diet can affect diversity at even the finest levels of taxonomy (De Filippis et al., 2019).  Interestingly enough, it also appeared to suggest how diets between different types of global populations, and in turn, the composition of gut microbiota, can further affect the ability of biological processes to consume certain types of food, most notably in non-Western countries like Italy (De Filippis et al., 2019).  Today, this species is now referred more commonly as *Segatella copri*, in line with the reclassification of many other species of this genus, from *Prevotella* to *Segatella* (Blanco-Miguez et al., 2023, Wang et al., 2026).

Studies such as these can prove most useful for explaining the essence of these connections between microbial taxonomy and lifestyle, in human health/disease as a whole and across populations (De Filippis et al., 2019).  However, it is also important that the composition and abundance of other taxa within the diversity should be investigated more thoroughly, especially at the species level, to understand exactly how they collectively support the necessary biological functions for food consumption and metabolism.  Hence, this taxonomic analysis seeks to expand on the significance of these connections, by applying different classification and diversity pipelines to compare how the abundances and diversity of certain microbiotic taxa in the human gut metagenomes differ across individuals with varying dietary patterns.  It will also attempt to explain how the presence and/or dominance of certain species across diets support digestive and metabolic functions and processes most associated with them.


## Methods
Raw short-read data for human gut microbiomes was randomly selected from samples used in a shotgun metagenomics study, with Illumina sequencing, from 74 Italians based on dietary habits, found on the National Center for Biotechnology Information (NCBI) database (University of Naples Federico II, 2017).  A total of six sets of paired-end samples were obtained, with three sets belonging to individuals with vegan diets and the other three belonging to individuals with omnivorous diets (University of Naples Federico II, 2017).  Quality control was conducted on each individual paired-end sample using FastQC (v0.12.1) on 12 threads per sample, to check for sources of contamination, which were then compiled into overall reports with MultiQC (v1.33) (Andrews, 2010, Ewels et al., 2016).  Adapter content was then trimmed off each paired-end sample using *fastp* on default parameters, before repeating quality control reports a second time and compressing them (.gz files) (Chen et al., 2018).  Taxonomic classification for each set of paired-end samples (.gz) was done using k-mer classification in Kraken2 (v2.1.6), using a confidence level of 0.15 to reduce the likelihood of false positives (Wood et al., 2019).  These computations were run remotely in the Narval Compute Canada cluster at Digital Research Alliance of Canada (DRAC) on 128G memory and 16 cores per task, with the most recent Kraken2 Standard database as of February 26, 2026 (Langmead, 2026).  Finally, taxonomic abundances for each sample were re-estimated at species level from classification genus reads, using Bayesian Re-estimation of Abundance with KrakEN (Bracken v3.0.1) at 8 threads per sample (Lu et al., 2017).  These Bracken species reports were subsequently compiled and converted into a BIOM file using *kraken-biom* (v1.2.0) for further taxonomic analysis (Dabdoub, 2016).

Continuing in R (v4.5.2), the BIOM report of all abundance tables from Bracken was imported into the working environment with the *phyloseq* package (v1.54.2) and the *biomformat* package (v1.38.3) (McMurdie & Holmes, 2013, McMurdie & Paulson, 2026).  Relative abundance of the top 20 most abundant taxa across all dietary groups and samples were compared and visualized with barplots using *phyloseq* (v1.54.2) and the *microbiome* package (v1.32.0), to explore the diversity of classifications and determine the most prominent microbial species in each set of metagenomes (McMurdie & Holmes, 2013, Lahti & Shetty, 2017).  Alpha and beta diversity measures were also both calculated using functions from the *phyloseq* package.  For alpha diversity, three different indices were calculated for each sample and visualized using dot plots.  Additionally, the mean &plusmn; standard Deviation was calculated overall for each of these indices in the omnivore and vegan diet groups.  The Shannon index was calculated to understand and visualize the ovarall taxa evenness in a sample, while the Simpson index was calculated to determine taxa dominance between groups and/or samples.  The Chao1 index was also calculated to evaluate the amount of taxa richness per sample.  Meanwhile, beta diversity was compared between vegan and omnivore diet groups using Jaccard and Bray-Curtis distance measures for dissimilarities in microbiomes.  Both measures of beta diversity were then visualized with PCoA plots.  Additionally, a Permutational Multivariate Analysis of Variance Based on Distance (PERMANOVA) statistical test was run (using the default 999 permutations) to test for significance in differences between the two groups using the *vegan* package (v2.7-3) (Oksanen et al., 2026).  Lastly, to determine the most significant differences in abundant taxa by diet, ANCOM-BC2 (v2.12.1) was used to conduct differential abundance analysis between the two groups (Lin & Peddada, 2020).  This analysis used the 5th percentile (0.05) of standard errors to define the *s0* constant, along with the default Holm-Bonferroni correction for multiple testing to reduce the false positive rate without being too fully conserative (Vinogradova et al., 2016).  Since adjusted p-values (*q*) returned a value of 1.00 for all samples after correction, though, the raw p-values (p) of taxa were used to subset the top 20 differentially abundant taxa to be visualized by effect size plot.

## Results
### Relative Taxonomic Abundance
Relative abundance for every sample was measured at the species level, with the top 20 taxa being ranked by highest abundance (Figure 1). In omnivore diets, the relative abundance between species looked more distributed for all the top 20 taxa.  Two species, though, *Bifidobacterium adolescentis* and *Faecalibacterium prausnitzii*, seemed to show higher relative abundance in their microbiomes.  In vegan diets, however, the abundance pattern seemed to be more distinct, as *Segatella copri*, the same species studied by De Filippis et al., showed a higher and more pronounced relative abundance across at least two of their samples (2019).

<img width="1902" height="1067" alt="Relative_Abundance" src="https://github.com/user-attachments/assets/f9a29ab6-2fc2-4d79-8464-fe120b41823c" />

**Figure 1**: Relative abundance of microbiome species from human gut metagenomes for Italian omnivore and vegan diets.  Classification and re-estimation steps in Kraken2/Bracken for taxa were done at the species level, and the top 20 taxa across all samples were selected and filtered for visualization based on highest abundance.

### Alpha Diversity
Alpha diversity was measured to quantify diversity within samples by examining the evenness, dominance, and richness of taxonomic species in the metagenomes.  Results were also compared between groups to evaluate how these measures differed in Italians based on omnivore and vegan diets.

#### Taxonomic Species Evenness - Shannon Index
The Shannon index was calculated to measure the amount of evenness in the overall distribution of all taxa for each sample (Figure 2).  When looking at the results, the mean &plusmn; SD evenness for omnivore diets was reported at about 3.67 &plusmn; 0.32, which was somewhat higher than the evenness for vegan diets, which was about 3.22 &plusmn; 0.18.  A closer look at Figure 2 showed that all three metagenomic samples belonging to the omnnivore diet expressed higher species evenness overall than the samples belonging to the vegan diet.  In the omnivorous group, the Shannon index ranged from about 3.45 to 4.04, a difference of 0.59.  This was due to one sample having a somewhat higher Shannon index as opposed to the other two samples, possibly highlighing a slightly skewed or inflated mean evenness.  In the vegan group, the Shannon Index ranged from 3.04 to 3.40, a difference of 0.36.  However, there was less difference observed between its three samples. 

<img width="1072" height="642" alt="Alpha_Diversity_Shannon_Index" src="https://github.com/user-attachments/assets/1ab3c12b-04e0-4808-ab19-0ba47197cd6f" />

**Figure 2**: Comparison of Shannon index for human gut metagenome samples between omnivore and vegan diets.  The mean &plusmn; SD evenness for omnivore diets was about 3.67 &plusmn; 0.32, whereas for vegan diets, it was about 3.22 &plusmn; 0.18.

#### Taxonomic Species Dominance - Simpson Index
Next, the Simpson index was used to assess dominance of certain taxa in abundance per sample, as well as to quantify the extent to which they did over taxa with average or lower abundance (Figure 3).  Due to this index being calculated using the *phyloseq* program, it represented the effect of dominance exhibited by one or more taxa.  Looking at all samples in both groups, the Simpson indices were all notably closer to 1 than to 0.  Additionally, they ranged only from about 0.83 to 0.97, a difference of 0.14, suggesting a higher overall dominant effect of one or more specific taxa, regardless of diet type.  Omnivorous diets had a mean &plusmn; SD index of 0.94 &plusmn; 0.03, while for vegan diets, it was about 0.89 &plusmn; 0.06, a smaller difference between the groups.  In Figure 3, there appeared to be less spread in the samples for the omnivore group, while the vegan group seemed to show slightly more consistent spread in its samples.  The Simpson indices for the omnivore group ranged from approximately 0.95 to 0.97, a very minimal difference of 0.02.  On the other hand, the vegan group had a range of about 0.83 to 0.94, a slightly greater difference of 0.11.

<img width="1055" height="597" alt="Alpha_Diversity_Simpson_Index" src="https://github.com/user-attachments/assets/3661af22-a2a8-4ef5-b8ca-a619eb7a9077" />

**Figure 3**: Comparison of Simpson index for human gut metagenome samples between omnivore and vegan diets.  The mean &plusmn; SD dominance for omnivore diets was about 0.94 &plusmn; 0.03, and for vegan diets, it was about 0.89 &plusmn; 0.06.

#### Taxonomic Species Richness - Chao1 Index
Finally, the Chao1 index helped to estimate and better establish overall species richness, which included both observed taxa and potentially unobserved taxa within each metagenomic sample (Figure 4).  Interestingly, unlike the first two indices, however, the vegan group actually had a higher mean and lower SD of 463 &plusmn; 68, as opposed to the omnivore group, which had a lower mean and higher SD at 441 &plusmn; 109.  Figure 4 also showed that, for the omnivorous diets, the sample which had a higher Shannon index than the other two samples similarly also had a higher Chao1 index than them, and the spread between it and the other samples was visually a lot greater this time.  Once again, this might have pointed towards a skewed effect on mean richness.  The range for this group was from about 341 to 557, a major difference of 216.  On the other hand, the vegan group also had one sample that was different from the other two, having a lower Chao1 index than them.  However, this group's range of indices was from 386 to 514, a lower difference of 128, possibly suggesting less spread or skewedness.

<img width="1037" height="552" alt="Alpha_Diversity_Chao1_Index" src="https://github.com/user-attachments/assets/0a2c1698-734d-4f8e-9095-1c3b08836d7e" />

**Figure 4**: Comparison of Chao1 index for human gut metagenome samples between omnivore and vegan diets.  The mean &plusmn; SD richness for omnivore diets was about 441 &plusmn; 109, while vegan diets had a higher value of 463 &plusmn; 68.

### Beta Diversity

<img width="1012" height="752" alt="Beta_Diversity_Jaccard" src="https://github.com/user-attachments/assets/6e931643-467f-492b-80cf-3a4bac85512f" />

<img width="995" height="707" alt="Beta_Diversity_Bray-Curtis" src="https://github.com/user-attachments/assets/9f8ae070-9503-416b-955e-363e4c61300a" />


### Differential Abundance - ANCOM-BC2
*Akkermansia muciniphila* & *Catenibacterium mitsuokai*

<img width="1920" height="1112" alt="Differential_Abundance_ANCOM-BC2" src="https://github.com/user-attachments/assets/523699aa-85eb-487d-a3e4-80e9059dc2b6" />



## Discussion
This taxonomic analysis of human gut metagenomics sought to provide a more insightful context of how the abundance and diversity of many different microbiota were connected to essential biological functions of vegan and omnivorous dietary habits.  However, considering that a low sample size was often a highly major limitation in the analyses, the scope of interpretation for these results, in turn, was often difficult to justify with full confidence, making them rather less reliably effective on their own.  Nonetheless, there were still some key findings that emerged which helped to better illustrate the relevance of these other abundant taxa to human health/disease.

Then, taking a further look into the alpha diversity of samples, it was determined that, on average, omnivorous samples appeared to have more species evenness in terms of abundance distribution, as well as strong dominance of one or more taxa, according to the Shannon and Simpson indices respectively.  However, it was the vegan samples that seemed to show more species richness on average when investigating the Chao1 index.

For beta diversity, 

Lastly, the differential abundance analysis revealed a mostly distinct separation in taxa abundance between the vegan and omnivorous groups.  Once again, though, they, along with all other taxa, were not proven statistically significant, so it was not possible to evaluate true significance of these differences.  All the same, though, it was still possible to understand how these findings were tied to existing trends or connections to the scientific literature on health and disease.  For vegan diets, the results suggested that *Akkermansia muciniphila* was the most enriched taxon for its group.  *A. muciniphila* is a very essential bacterium that resides in the gastrointestinal gut, responsible for the degradation of a protective mucus layer known as mucin (Collado et al., 2007).  It has been found to inhabit and colonize this area of the digestive system as early as infant stages, and its role in mucin degradation helps to release nutrients from the layer that are essential for the growth of other necessary microbiota in this environment (Collado et al., 2007).  Mucin itself, .  On the other hand, for omnivorous diets, the results suggested *Catenibacterium mitsuokai* as the most abundant in comparison.

## References
Andrews, S. (2010). FastQC. A quality control tool for high throughput sequence data. *Babraham Bioinformatics*. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Blanco-Miguez, A., Galvez, E. J. C., Pasolli, E., De Filippis, F., Amend, L., Huang, K. D., Manghi, P., Lesker, T. R., Riedel, T., Cova, L., Puncochar, M., Thomas, A. M., Valles-Colomer, M., Schober, I., Hitch, T. C. A., Clavel, T., Berry, S. E., Davies, R., Wolf, J., Spector, T. D., Overmann, J., Tett, A., Ercolini, D., Segata, N., & Strowig, T. (2023). Extension of the Segatella copri complex to 13 species with distinct large extrachromosomal elements and associations with host conditions. *Cell Host Microbe*, *31*(11), 1804-1819.e9. https://doi.org/10.1016/j.chom.2023.09.013

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, *34*(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

Collado, M. C., Derrien, M., Isolauri, E., de Vos, W. M., & Salminen, S. (2007). Intestinal Integrity and Akkermansia muciniphila, a Mucin-Degrading Member of the Intestinal Microbiota Present in Infants, Adults, and the Elderly. *Applied and Environmental Microbiology*, *73*(23), 7767-7770. https://doi.org/10.1128/AEM.01477-07

Dabdoub, S. M. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2). https://github.com/smdabdoub/kraken-biom

De Filippis, F., Pasolli, E., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, *25*(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, *32*(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

Lahti, L., Shetty, S., et al. (2017). Tools for microbiome analysis in R. Microbiome package version 1.32.0. http://microbiome.github.com/microbiome 

Langmead, B. (2026). Index zone - Kraken 2 / Bracken Refseq indexes. https://benlangmead.github.io/aws-indexes/k2

Lin, H. & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*, *11*, 3514. https://doi.org/10.1038/s41467-020-17041-7

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, *3*, e104. https://doi.org/10.7717/peerj-cs.104

McMurdie, P. J. & Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLOS One*, *8*(4), e61217. https://doi.org/10.1371/journal.pone.0061217

McMurdie P, Paulson J (2026). biomformat: An interface package for the BIOM file format [R package version 1.38.3]. https://github.com/joey711/biomformat/

Oksanen, J., Simpson, G., Blanchet, F., Kindt, R., Legendre, P., Minchin, P., O'Hara, R., Solymos, P., Stevens, M., Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B., Borcard, D., Borman, T., Carvalho, G., Chirico, M., De Caceres, M., Durand, S., Evangelista, H., FitzJohn, R., Friendly, M., Furneaux, B., Hannigan, G., Hill, M., Lahti, L., Martino, C., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., Smith, T., Stier, A., Ter Braak, C., & Weedon, J. (2026). vegan: Community Ecology Package [R package version 2.7-3]. https://doi.org/10.32614/CRAN.package.vegan
  
University of Naples Federico II. (2017). SRX4967482: gut metagenome of vegan subject from Bari (Italy), subject ID 01BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967482[accn]

University of Naples Federico II. (2017). SRX4967479: gut metagenome of vegan subject from Parma (Italy), subject ID 26PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967479[accn]

University of Naples Federico II. (2017). SRX4967478: gut metagenome of omnivore subject from Bari (Italy), subject ID 03BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967478[accn]

University of Naples Federico II. (2017). SRX4967471: gut metagenome of vegan subject from Turin (Italy), subject ID VOV56. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967471[accn]

University of Naples Federico II. (2017). SRX4967459: gut metagenome of omnivore subject from Turin (Italy), subject ID VOV77. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967459[accn]

University of Naples Federico II. (2017). SRX3463060: gut metagenome of omnivore subject from Parma (Italy), subject ID 37PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX3463060[accn]

Vinogradova, E., Kushugulova, A., Kozhakhmetov, S., & Baltin, M. (2026). The Use of Confidence Intervals in Differential Abundance Analysis of Microbiome Data. *Applied Microbiology*, *6*(1), 7. https://doi.org/10.3390/applmicrobiol6010007

Wang, S., Zhou, T., Wang, X., Zhao, J., & Wang, X. (2026). Bridging the gap: Prevotella/Segatella's impact on gut barrier function and advanced cultivation strategies to realize the uses in gut health. *Gut Microbes*, *18*(1), 2638001. https://doi.org/10.1080/19490976.2026.2638001

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, *20*, 257. https://doi.org/10.1186/s13059-019-1891-0
