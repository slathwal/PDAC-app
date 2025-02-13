## Introduction
This application contains analysis of Pancreatic Ductal Adenocarcinoma 
(PDAC) samples from the Cancer Genome Atlas (TCGA). A brief description of the analyses included in the application is provided below. For more information about PDAC data from TCGA, please visit the PDAC official publication by TCGA, [Integrated Genomics Characterization of Pancreatic Ductal Adenocarcinoma](https://www.cell.com/cancer-cell/fulltext/S1535-6108(17)30299-4?).

## Demographic Overview
This page presents summary statistics of key demographic variables in the study - age, gender, and race

## Censoring Plot
This page shows the censor status of each subject, broken down by gender. A subject is said to be censored at a particular time if they do not show the event of interest by a given time and are not followed further in the study.

## Kaplan Meier plot for survival
This page shows the Kaplan Meier survival curves for subjects in the study. The page includes a dropdown that can be used to select variables. When a variable is selected, the survival curves are split by the values of the selected variable and a table containing p-values for pariwise log-rank test between all the values of the selected variable is shown.

## Principal Component Analysis
This page shows the first and second principal components of the proteomics data on the subjects. The page contains a dropdown that can be used to select variables to colour the subjects by available metadata.

## Tumor Purity Analysis
This page visualizes the relationships between estimated purity of the tumor samples and molecular subtype assigned to samples based on previous gene expression research by [Moffitt et al.](https://doi.org/10.1038/ng.3398), [Bailey et al.](https://doi.org/10.1038/nature16965), and [Collisson et al.](https://doi.org/10.1038/nm.2344)

## Acknowledgement
The design of the application is inspired by [Pilot 2](https://rconsortium.github.io/submissions-wg/pilot2.html) of the R Submissions Working Group at R Consortium.