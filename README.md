# Deep-Learning-in-rice-for-prediction

Using Bayesian linear models, we have shown that Transposable Insertion Polymorphisms (TIPs) can improve prediction ability in genomic prediction of complex agronomic traits in rice over standard approaches based exclusively on Single Nucleotide Polymorphisms (SNPs). However, TIPs are not the only structural variation in the genome. Structural variations (SVs) such as deletions, inversions, and duplications are prevalent in the genome and they play an important role in plant evolution. Here, we investigate whether combining the structural and nucleotide genome-wide variation can improve the prediction ability of traits when compared to using only SNPs. For the purposes of the study, four important agronomic traits were used from 738 rice accessions in total, originating from five different rice population groups (Aus, Indica, Aromatic, Japonica, and Admixed). We assess prediction accuracy by applying cross-validation under two different strategies. In the first strategy, we used a k-fold cross-validation producing ten partitions from the whole population. In the second strategy, we followed an across-population scenario that predicted Aromatic and Admixed accessions from the rest of the populations. In each scenario, the performance of BayesC and a Bayesian Reproducible Kernel Hilbert space regressions are compared to Deep Learning networks. We investigated the prediction ability of DL using two different widely used architectures, a Multilayer Perceptron, and a Convolution Neural Network. Then we further explored their performance by using various marker input strategies. We found that exploiting structural and nucleotide variation improves the prediction ability of complex traits in rice. Also, our results suggested that DL models outperform in 75% of the studied cases. Finally, DL constantly improved the prediction ability of binary traits against the Bayesian models.


## FILES
### SCRIPTS
* BAYESC_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using BayesC from BGLR package.
* RKHS_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using RKHS from BGLR package.
* Hypermodel_CNN.py: Python script for implementing Convolutional Neural Network generated for genomic prediction.
* Hypermodel_MLP.py: Python script for implementing Multilayer Perceptron Neural Network generated for genomic prediction.
* Hypermodel_Multiple_inputs.py: Python script for implementing Multilayer Perceptron Neural Network using multiple inputs for genomic prediction.
* ngsLD_script.sh: Bash script for running ngsLD software.

### RESULTS_TABLES
* accuracy.dl.results: Table of accuracy values resulted by Deep Learning analysis.
* accuracy.linear.results: Table of accuracy values results by Bayesian linear models.
* corr.dl.results: Table of correlation values resulted by Deep Learning analysis.
* corr.linear.results: Table of correlation values results by Bayesian linear models.
* loss.dl.results: Table of loss values resulted by Deep Learning analysis.
* loss.linear.results: Table of loss values results by Bayesian linear models.
* unique.linked.snps.across.all.the.analysis: TXT file with the name and positions of all the linked snps found in the analysis.

### DATA
* final_deletions.RData: Genotypes of deletions in R data frame format.
* final_duplications.RData: Genotypes of duplications in R data frame format.
* final_inversions.RData: Genotypes of inversions in R data frame format.
* final_mitedtx.RData: Genotypes of mitedtx in R data frame format.
* final_rlxrix.RData: Genotypes of rlxrix in R data frame format.
* final_snps.RData: Genotypes of snps in R data frame format.
