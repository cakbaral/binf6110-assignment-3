# binf6110-assignment-3

## Introduction

## Methods
Raw short-read data for human gut microbiomes was randomly selected from samples used in a shotgun metagenomics study on 74 Italians based on dietary habits, on the National Center for Biotechnology Information (NCBI) database (University of Naples Federico II, 2017).  A total of six sets of paired-end samples were obtained, with three sets belonging to individuals with vegan diets and the other three belonging to individuals with omnivorous diets (University of Naples Federico II, 2017).  Quality control was conducted on each individual paired-end sample using FastQC (v0.12.1) to check for sources of contamination, and the reports were compiled into overall reports with MultiQC (v1.33) (, Ewels et al., 2016).  Adapters were then automatically detected and trimmed off each paired-end sample using *fastp* before repeating quality control reports a second time ().  For taxonomic classification of samples, each set of paired-end samples was classified through Kraken2 (v2.17.1) using the Standard-16 database ().  The using Bayesian Reestimation of Abundance with Kraken (Bracken) (). 

Alpha and beta diversity measures were both calculated using the *vegan* package .  For alpha diversity, the Shannon and Simpson index were used to .  Meanwhile,  

## Results

## Discussion

## References
Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, *32*(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

University of Naples Federico II. (2017). SRX4967482: gut metagenome of vegan subject from Bari (Italy), subject ID 01BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967482[accn]

University of Naples Federico II. (2017). SRX4967479: gut metagenome of vegan subject from Parma (Italy), subject ID 26PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967479[accn]

University of Naples Federico II. (2017). SRX4967478: gut metagenome of omnivore subject from Bari (Italy), subject ID 03BA. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967478[accn]

University of Naples Federico II. (2017). SRX4967471: gut metagenome of vegan subject from Turin (Italy), subject ID VOV56. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967471[accn]

University of Naples Federico II. (2017). SRX4967459: gut metagenome of omnivore subject from Turin (Italy), subject ID VOV77. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX4967459[accn]

University of Naples Federico II. (2017). SRX3463060: gut metagenome of omnivore subject from Parma (Italy), subject ID 37PR. *NCBI*. https://www.ncbi.nlm.nih.gov/sra/SRX3463060[accn]

De Filippis, F., Pasolli, E., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, *25*(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004
