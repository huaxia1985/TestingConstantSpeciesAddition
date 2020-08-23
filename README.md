# TestingConstantSpeciesAddition
R codes to generate null models and testing model adequacy in phylogenetic analyses of temporal patterns in species addition
See details in paper "Modelling colonization rates over time: generating null models and testing model adequacy in phylogenetic analyses of species assemblages" by Xia Hua and Lindell Bromham
"code for colonization timing" is the executable code to generate Figure S1-S3 in the paper.
"code for tip branch length" is the executable code to generate Figure S4-S6 in the paper.

Folder 'Madagascar squamates' includes all the codes for testing constant colonisation rates in the Madagascar squamates case study.
main.R is the executable code.
colcal.R infers the probability distribution of colonisation times from assemblage phylogeny
colnllcal.R derives the null distribution of colonisation times
colsigtest.R tests whether the inferred colonisation times are random draws from the null distribution of colonisation times
extractcome.R extract subtree information from the assemblage phylogeny
likcal.R calculates the likelihood of the subtrees given a GeoSSE model 
