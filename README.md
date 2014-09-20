iPaD
====

Authors: 

Cong Li,      cong.li@yale.edu
Can Yang,     eeyang@hkbu.edu.hk
Greg Hather,  ghather@gmail.com
Ray Liu,      ray.liu@takeda.com
Hongyu Zhao,  hongyu.zhao@yale.edu

====

This is a package written in Matlab that performs drug-pathay association analysis on paired transcription/drug sensitivity profile data using an efficient bi-convex optimization algorithm. 

It requires the following four matrices as input: 1) An N by G transcription profile matrix Y1; 2) An N by D drug sensitivity matrix Y2; 3) A P by G pathway-gene relationship indicator matrix L1; 4) A P by D drug-pathway association indicator matrix L2_prior. N is the number of cell lines (samples), G is the number of genes, D is the number of drugs and P is the number of pathways. Of course, the samples, genes, drugs and pathways have to be in the same order across these matrices. 

====

This repository contains the following files:

iPaD.m        This is the main function of the iPaD package.
iPaD_cv.m     Perform cross-validation to choose an appropriate penalty parameter.
iPaD_permu.m  Perform permutaton test to calculate the p-values for the drug-pathway associations.

(initialize_X.m LassoSolver.m update_B1.m update_B2.m update_X.m) 
These files are all sub-routines of the iPaD main function.

quality_control.m   Perform the two following quality control steps for the input data set: 1) remove genes or drugs that have less than three unique values; 2) merge pathways that have identical member genes.

CCEL_analysis.m   Code for analyzing an example real data set - the Cancer Cell Line Encyclopedia (CCLE) data set. Can be used as a vignette for the usage of the iPaD package.
(CCLE_L1.txt CCLE_L2_prior.txt CCLE_L2_validate.txt CCLE_Y1.txt CCLE_Y2.txt CCLE_genes.txt CCLE_drugs.txt CCLE_pathways.txt)
These files are the CCEL data set.

NCI60_analysis.m   Code for analyzing another example real data set - the NCI-60 data set. Also can be used as a vignette for the usage of the iPaD package.
(NCI60_L1.txt NCI60_L2_prior.txt NCI60_L2_validate.txt NCI60_Y1.txt NCI60_Y2.txt NCI60_genes.txt NCI60_drugs.txt NCI60_pathways.txt)
These files are the NCI-60 data set.

