# HNCRM
Hierarchical Normalized Completely Random Measures to Cluster Grouped Data

# Authors

>  Raffaele Argiento, ESOMAS Department, University of Torino and Collegio Carlo Alberto, Torino, Italy;

> Andrea Cremaschi, Department of Cancer Immunology, Institute of Cancer Research, Oslo University Hospital, Oslo, Norway
and Oslo Centre for Biostatistics and Epidemiology, University of Oslo, Oslo, Norway;

> Marina Vannucci, Department of Statistics, Rice University, Houston, TX, USA.

# Description
A Bayesian nonparametric model for clustering grouped data is described. We adopt a hierarchical approach: at the highest level, each group of data is modelled according to a mixture, where the mixing distributions are conditionally independent normalized completely random measures (NormCRMs) centered on the same base measure, which is itself a NormCRM. The discreteness of the shared base measure implies that the processes at the data level share the same atoms. This desired feature allows to cluster together observations of different groups. We obtain a representation of the hierarchical clustering model by marginalizing with respect to the infinite dimensional NormCRMs. We investigate the properties of the clustering structure induced by the proposed model and provide theoretical results concerning the distribution of the number of clusters, within and between groups. Furthermore, we offer an interpretation in terms of generalized Chinese restaurant franchise process, which allows for posterior inference under both conjugate and non-conjugate models. We develop algorithms for fully Bayesian inference and assess performances by means of a simulation study and a real-data illustration.

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

