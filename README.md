# Evaluation of Deep Learning for predicting rice traits using structural and single-nucleotide genomic variants
Ioanna-Theoni Vourlaki1, Sebastián E. Ramos-Onsins1,Miguel Pérez-Enciso1,2,3* Raúl Castanera1*

1Centre for Research in Agricultural Genomics CSIC-IRTA-UAB-UB, Campus UAB, Edifici CRAG, Bellaterra, Barcelo-na 08193, Spain,
2Catalan Institute for Research and Advanced Studies (ICREA), Barcelona, Spain,
3Universitat Autónoma de Barcelona, Barcelona, 08193, Spain.
For correspondence: mperezenciso@gmail.com, raul.castanera@cragenomica.es
 
Structural variants (SVs) such as deletions, inversions, duplications, and Transposable Element (TE) Insertion Polymorphisms (TIPs) are prevalent in plant genomes and have played an important role in evolution and domestication, as they constitute a significant source of ge-nomic and phenotypic variability. Nevertheless, most methods in quantitative genetics that focus on plant crop improvement, such as genomic prediction, consider Single Nucleotide Pol-ymorphisms (SNPs) as the only type of genetic marker. Here, we used rice to investigate whether combining the structural and nucleotide genome-wide variation can improve prediction ability of traits when compared to using only SNPs. Moreover, we also examine the potential ad-vantage of Deep Learning (DL) networks over Bayesian Linear models, which have been widely applied in genomic prediction. Specifically, the performance of BayesC and a Bayesian Repro-ducible Kernel Hilbert space regressions were compared to two different DL architectures, the Multilayer Perceptron, and the Convolution Neural Network. We further explore their predic-tion ability by using various marker input strategies and found that exploiting structural and nucleotide variation improves prediction ability on complex traits in rice. Also, DL models out-performed Bayesian models in 75% of the studied cases. Finally, DL systematically improved prediction ability of binary traits against the Bayesian models.


## FILES
### SCRIPTS
* BAYESC_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using BayesC from BGLR package.
* RKHS_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using RKHS from BGLR package.
* Hypermodel_CNN.py: Python script for implementing Convolutional Neural Network generated for genomic prediction.
* Hypermodel_MLP.py: Python script for implementing Multilayer Perceptron Neural Network generated for genomic prediction.
* Hypermodel_Multiple_inputs.py: Python script for implementing Multilayer Perceptron Neural Network using multiple inputs for genomic prediction.
* ngsLD_script.sh: Bash script for running ngsLD software for each trait and partition.

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
