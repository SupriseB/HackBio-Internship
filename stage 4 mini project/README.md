<!--StartFragment-->

**Machine Learning modeling and Biomarker Discovery in Gliomas Based on Methylation Data: Insights into IDH Status Classification**



****Contributors:**
**
@Oluwakamiye, @Sapi, @Mojoyin, @Berkay, @Ifeoluwa01.

Code link: <https://github.com/hberkayc/HackBio_Stage4>

**Introduction**

Glioma is a type of brain tumor which affects glial cells. Isocitrate dehydrogenase (IDH) gene, can be used as a biomarker for the classification of gliomas since it is closely associated with tumor development and prognosis (Agnihotri et al., 2014).

The aim of this task was to reproduce the unsupervised clustering of gliomas based on methylation levels to classify six distinct IDH statuses and to identify expression biomarkers that can differentiate between pairs of IDH classes, providing deeper insights into the genetic and molecular underpinnings of glioma subtypes**.**

 ****

**Methods**

The dataset for this task was obtained from the TCGA database using TCGABiolinks in R. TCGA DNA methylation data for 516 glioma samples were preprocessed, normalized, and filtered. Differential methylation analysis was performed using limma to compare wildtype and mutant IDH statuses and functional enrichment analysis using GO and KEGG. This was followed by feature selection, the application of KNN and Random-Forest and K-means for clustering, followed by an evaluation of the model performance

We performed three machine learning algorithms: KNN, Random Forest and K-means clustering- Even though the question says we should run KNN algorithm, we decided to run Random forest and also K-means clustering for model comparisons. KNN and Random forest gave an accuracy of 100% and K-means clustering gave an accuracy of 99.81%

To ensure robustness and avoid overfitting, we employed 5-fold cross-validation. In this process, the dataset was randomly divided into five subsets (folds), with four subsets used for training and one subset for validation in each iteration. This process was repeated five times so that each fold served as a validation set exactly once.

\
\
\
\


**Results**

**Biomarker analysis**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdib8rgDoEcanVnjiI9ZdNg5sLg5-qBN-JISaDxjPQ4Tk7gyF3R0jgTVnUY9Ow9DtrX880jqjq8iHI5Xc1GUuwCRBUGOIaTB6xHpW8AK08X2zTQL3xrqKXuqcZNwUt51-caEbw8kwNrwklnmWLJOKSMoiG8?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.1: Post-normalization Boxplot**

To normalize the “Methylation Beta Values”, we used the “quantile” method. 

\


****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeuZGGWWIi5TwwVoePFreq2zbZCmTnHBn6AndiFI9o9q8hrPldIdKAs4nHvMfnV6-bnrCj_grElDtmBUbkCV5EXyT9VO9UjGB8jXJ3MjSpvmfhmItUTW3zufWlpxkDgH_M78XMCbR-mMDefRWRS7IM52ZKd?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.2: Heatmap of the probes (rows) and samples (columns)**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXd8dku22D9idaa5BibkAKWbYa57cds6f-T2_OZI9s1uwa_WVshMYVAGuOFwUm8D644JM63d_5i6G-UeIPqySeWkLVnngC_YJofxLPkIoDVUWuZUp9Poehp6nj-o7KIRaorUsMksJjXLLvBeob4hEzWVKFvG?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.3: Heatmap of the top 50 probes**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc6A9DjxcXfJRRft8Z1IscemQk412BLLHXX4O2g32Zeqma2LKqJ1hHbB7LwPhVi_Z9WK6qJrDo4Q5SjLfh8ONzWK0KwOSxCFKtwb1EavOhQ6P7utHaYIh3_zv8ygouIDtRb37lZAwXHSBtbplGvWc4SMPwl?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.4: Volcano plot of the differential methylation analysis**

To detect the significant probes, we set the logFC\_cutoff <- 0.4, pvalue\_cutoff <- 0.05.  

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfNE8y7F0h7MUlUpk6L_UzYKGIIO2GPEyFIwQKX8VxfXL1d1t64cF5Y4gM1A-D6A6mKVIXuUmRkpDOlKvgBD1jxqWTKXRs63yeLHyauiMnpfG-W5goQTnijUIZ47Jiv1ZrXfI-ufuCZEiTaXmIxaZEvS8Q?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.5: Dotplot of the Gene Ontology (GO) enrichment analysis of all probes (corresponding to their genes)**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcaDVJ14G1EIB7OXbgwDkL_4pzQxwH56NyFmjy25xQ_4sf76ix7klqVEjG6U4HcblJmlmTkAfKzTFBquPwVgGjeT108Aaje2GAfFCudwxF5x-V-x3ttsWlRNwv7RHzP6xduvavIMI8Iaq_Lw3c6U-s2PDc?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.6: cnetplot of the Gene Ontology (GO) enrichment analysis of all probes (corresponding to their genes)**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXds92LCgY1rvVDeh6ahEFdsiHFQTNd_IZJy-ZqKyAmsSM9y5DNpGzen3UU7X2lJqnl63Rko3SuUNWInQ5Van8tL2yJ28Zq84O087VHOFmruCSJQ8xWtSCRWcQdraStx_Em0LfNuwLJIWIpYg8YOd1WlfNrc?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.7: Dotplot of the Gene Ontology (GO) enrichment analysis of top 50 probes (corresponding to their genes)**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdkMPwr8iB6qpal8ihc2q0q_GKX0Nqb_0DFjwlXqSLvefsxVyjs-EAGM4hXvMtnyd_ZslkKaWXXMO53vxzdF-EpZYVYRkppgKxnwpcd3O-ScLdNBP-6bXx860iA6ePZcpk2uF18l0KCI2cO-_Ck8JUCHr1q?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.8: cnetplot of the Gene Ontology (GO) enrichment analysis of top 50 probes (corresponding to their genes)**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdO6SjFQ_tVUHxdGnD_g2rsRWWvA-LoU4Mq-2AthHu_uKRBIYc6ZkVA6iEgz-s4BQw1ro8DBIFgThnmEXCmYr-HDzoFY7gSj9lymlTV4yDOFuE3BUpocf1ZhvTOOuttyVdwYOA7ZjxuTknM9sXb12Lb0jsK?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure.9: Dotplot of the KEGG enrichment analysis of top 50 probes (corresponding to their genes)**

For the KEGG analysis, we couldn’t get a result for all of the probes. 

 **Machine Learning Results**

Random Forest ROC AUC scores for each fold: \[1. 1. 1. 1. 1.]

Random Forest Mean ROC AUC: 1.0000

KNN ROC AUC scores for each fold: \[1.     0.97771739 1.     1.     1.    ]

KNN Mean ROC AUC: 0.9955

 

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeQlVTEughlb-0zTowkHvJqTEbl0Zz9D7VabSXscae2tkUboEa6Ymdwu_B2M2OD5uW-ctbJg6YeS8PJXde-K38E7D_RNgBdUjoJGAcjZXxZPpVt20OWnFU4E7qbNmDc16-fbKUmon966i7jbC77oUqSEHqI?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure 10:  KNN Confusion Matrix**

KNN Precision: 1.0

KNN Recall: 1.0

Random Forest Model Accuracy: 100.00%

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcmRzq6oU1GCB9KfAf7T5xRH0xUJWSEU_lV0fgSuS5Gxmh44zFOb-do3kWav5RCExJ1xSZlsymyrv-766ZIw4Epau2ha5VXhOV6rsdvyhBt0_oO_ScjSrEF-NlI2v-Rr2wlYp-s6FQ4lH5Pxgzh4mcNMCZ3?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure 11: Random forest confusion matrix**

Random Forest Precision: 1.0

Random Forest Recall: 1.0

Random Forest Model Accuracy: 100.00%

\


****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfK37J7sL5NtcPqVCG1sdQwV3a0lG3PQW-vvcWBhMZXhMmA_wGaLFoXWgjQNaew-_ads3H4ZA8DCs9OpUmtN1Jtk4EeKBt35CiMkZvgDJI57x6tt_5J6sNMZT0Zi3JZuD51R2tvdfAGOBzQZ0jvdM_gOjZm?key=HKP7dWO6HoopXIJggI4OJA)****

**Figure 12: Elbow method for evaluation of optimal  k-means**

The Elbow method helps us to determine the optimal number of clusters by plotting the inertia (or WCSS) for different numbers of clusters. When the reduction in inertia starts to slow down (forming an "elbow"), that's often the optimal number of clusters.

With our known class labels (IDH\_status), we compared the clusters to the known labels to see how well they corresponded. We calculated metrics like the Adjusted Rand Index (ARI) or Normalized Mutual Information (NMI) to compare the clustering results to the ground truth.

ARI compares the clustering result with ground truth labels by evaluating the number of pairs of points that are either correctly clustered together or correctly separated. It adjusts for chance by considering how well random cluster assignments would perform.

Range: ARI scores range from -1 to 1.

1:Perfect match between the predicted clusters and true labels (all points are clustered perfectly).

0: Clustering result is random and doesn't provide any better-than-chance clustering.

Negative values: Worse than random clustering.

What a high ARI score indicates: A high ARI score means that the clustering result is very similar to the ground truth labels,with minimal random or incorrect assignments. It shows that the model correctly grouped similar points together and separated dissimilar ones.

NMI is based on the concept of information theory. It measures the amount of information shared between the predicted clusters and true labels. The higher the shared information, the better the clustering matches the true labels.

Range: NMI scores range from 0 to 1. 1: Perfect correlation between the predicted clusters and true labels (the clusters provide all the information needed to know the true labels). 0: No mutual information between the predicted clusters and true labels (the clustering result doesn't give any useful information about the true labels). What a high NMI score indicates: A high NMI score indicates that the predicted clusters align well with the true labels, with minimal entropy or uncertainty. The clusters are informative and reflect the true class distribution well.

**Conclusions**

Differential methylation analysis revealed significant differential methylation between mutant and wild-type (WT) groups. Mutants showed widespread hypermethylation and hypomethylation compared to WT. These findings highlight potential biomarkers or pathways underlying group differences.

The most enriched and statistically significant biological processes in this dataset are related to cellular responses to inorganic substances and protein localization, which might be important in understanding the underlying biological context.

The machine learning analysis showed strong performance of KNN, Random Forest, and K-means models, with Random Forest being most reliable.

\
\
\
\
\
\
\
\
\


**References**

- AGNIHOTRI, S., ALDAPE, K. D. & ZADEH, G. 2014. Isocitrate dehydrogenase status and molecular subclasses of glioma and glioblastoma. _Neurosurgical focus,_ 37, E13.

- Carlson M. (2024). \_org.Hs.eg.db: Genome-wide annotation for human\_. R package version 3.19.1.

- Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2024). “edgeR 4.0: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” _bioRxiv_.[ doi:10.1101/2024.01.21.576131](https://doi.org/10.1101/2024.01.21.576131).

- Gao, Y., Zhang, H. & Tian, X. (2024). Integrated analysis of TCGA data identifies endoplasmic reticulum stress-related lncRNA signature in stomach adenocarcinoma. _Oncologie_, _26_(2), 221-237.[ https://doi.org/10.1515/oncologie-2023-0394](https://doi.org/10.1515/oncologie-2023-0394)

- Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters.OMICS: A Journal of Integrative Biology 2012, 16(5):284-287

- McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” _Nucleic Acids Research_, 40(10), 4288-4297.[ doi:10.1093/nar/gks042](https://doi.org/10.1093/nar/gks042)

- Mounir, Mohamed, Lucchetta, Marta, Silva,   T, Olsen, Catharina, Bontempi, Gianluca, Chen, Xi, Noushmehr, Houtan, Colaprico, Antonio, and Papaleo, Elena (2019). “New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.” _PLoS computational biology_, 15(3), e1006701.

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47.[ doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007). 

- Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a bioconductor package for differential expression analysis of digital gene expression data.” _Bioinformatics_, 26(1), 139-140.[ doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).

- S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, doi:10.1038/s41596-024-01020-z

- Silva, C. T., Colaprico, Antonio, Olsen, Catharina, D'Angelo, Fulvio, Bontempi, Gianluca, Ceccarelli, Michele, Noushmehr, Houtan (2016). “TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages.” _F1000Research_, 5.

- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). \_dplyr: A Grammar of Data Manipulation\_. R package version 1.1.4, <[https://CRAN.R-project.org/package=dplyr](https://cran.r-project.org/package=dplyr)>.

<!--EndFragment-->
