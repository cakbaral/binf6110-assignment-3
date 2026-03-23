# binf6110-assignment-3

## Introduction

## Methods
Raw short-read data for human gut microbiomes was randomly selected from samples used in a shotgun metagenomics study on 74 Italians based on dietary habits, on the National Center for Biotechnology Information (NCBI) database (University of Naples Federico II, 2017).  A total of six sets of paired-end samples were obtained, with three sets belonging to individuals with vegan diets and the other three belonging to individuals with omnivorous diets (University of Naples Federico II, 2017).  Quality control was conducted on each individual sample using FastQC and MultiQC, and adapters were trimmed with *fastp* ().  For taxonomic classification of samples, each set of paired-end samples was classified through Kraken2 (v2.17.1) using the Standard-16 database ().  The using Bayesian Reestimation of Abundance with Kraken (Bracken) (). 

Alpha and beta diversity measures were both calculated using the *vegan* package .  For alpha diversity, the Shannon and Simpson index were used to .  Meanwhile,  

## Results

## Discussion

## References
University of Naples Federico II. (2017). gut metagenome of vegan subject from Parma (Italy), subject ID 26PR. *NCBI*. [https://www.ncbi.nlm.nih.gov/sra/?term=SRP126540](https://www.ncbi.nlm.nih.gov/sra/SRX4967479[accn])
De Filippis, F., Pasolli, E., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, *25*(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004
