# Evaluation of Deep Learning for predicting rice traits using structural and single-nucleotide genomic variants
Ioanna-Theoni Vourlaki1, Sebastián E. Ramos-Onsins1,Miguel Pérez-Enciso1,2,3* Raúl Castanera1*

1Centre for Research in Agricultural Genomics CSIC-IRTA-UAB-UB, Campus UAB, Edifici CRAG, Bellaterra, Barcelona 08193, Spain,
2Catalan Institute for Research and Advanced Studies (ICREA), Barcelona, Spain,
3Universitat Autónoma de Barcelona, Barcelona, 08193, Spain.
For correspondence: mperezenciso@gmail.com, raul.castanera@cragenomica.es
 
Abstract
Background
Structural genomic variants (SVs) are prevalent in plant genomes and have played an important role in evolution and domestication, as they constitute a significant source of genomic and phenotypic variability. Nevertheless, most methods in quantitative genetics focusing on crop improvement, such as genomic prediction, consider only Single Nucleotide Polymorphisms (SNPs). Deep Learning (DL) is a promising strategy for genomic prediction, but its performance using SVs and SNPs as genetic markers remains unknown.
Results
We used rice to investigate whether combining SVs and SNPs can result in better trait predic-tion over SNPs alone and examine the potential advantage of Deep Learning (DL) networks over Bayesian Linear models. Specifically, the performances of BayesC and a Bayesian Repro-ducible Kernel Hilbert space (RKHS) regression, were compared to those of two different DL architectures, the Multilayer Perceptron, and the Convolution Neural Network. In the case of RKHS, models were implemented considering either additive effects exclusively or non-additive additionally.  We further explore their prediction ability by using various marker in-put strategies. We found that exploiting structural and nucleotide variation slightly improved prediction ability on complex traits in 87% of the cases, and that DL models outperformed Bayesian models in 75% of the studied cases. Finally, DL systematically improved prediction ability of binary traits against the Bayesian models.



## FILES
### SCRIPTS
* BAYESC_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using BayesC from BGLR package.
* RKHS_MODEL.R: R script for genomic prediction in 11 different cross-validation scenarios using RKHS from BGLR package.
* RKHS_MODEL_NON_ADDITIVE_EFFECTS.R: R script for genomic prediction in 11 different cross-validation scenarios using RKHS incorporating additive and non-additive effects.
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
