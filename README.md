# HNCRM
Hierarchical Normalized Completely Random Measures to Cluster Grouped Data

# Contents
1) N_Ind_NGG.m
2) N_HNGG.m
3) School_HNGG.m

Each of the files performs MCMC sampling according to the models described in the manuscript, for the simulated data (N_Ind_NGG.m and N_HNGG.m), or for the School data (School_HNGG.m).

The MCMC output is saved, and can be used to produce the figures presented in the paper.

The code producing the figures follows the main MCMC algorithm part.

#################
1) N_Ind_NGG.m
Matlab code for simulated data (d = 2) fitted using the conjugate Normal-inverse gamma model and independent NGG processes for each group.

#################
2) N_HNGG.m
Matlab code for simulated data (d = 2) fitted using the conjugate Normal-inverse gamma model and hierarchical NGG process.

#################
3) School_HNGG.m
Matlab code for the School data (Hoff, 2009) fitted using the a non-conjugate Normal-inverse gamma model and hierarchical NGG process.

