# DrugInteractionsPipeline
Source code for data analysis - Brochado et al Nature 2018

This folder contains source code implemented and used to analyze the data published by Brochado et al, Nature 2018 (DOI:).
RStudio Version 1.0.136, and R version 2.15.1 were used.
This analysis pipeline can in principle run without a high-performance computing system. However, the authors recommend access to a high-performance computing system to run specified parts of the pipeline.
Furthermore, the authors advice caution in re-using this pipeline for other datasets than the one it was built for, due to potential divergences on data quality.

## The pipeline
The pipeline is rooted at *DrugInteractionsPipeline* and split in 3 parts:
1) Analysis of drug interactions at the root
2) Sensitivity analysis
3) Interaction-scores calculator (next update)

They should be used following the aforementioned order, in order to ensure that all functions and data dependences are satisfied.
The authors recommend a careful look into the directory paths for input-output file allocation in order to ensure that all paths are properly set, before running the entire pipeline.
All figures related to data analysis relevant for the publication mentioned above were generated with these scripts.

## Input files & running instructions
The source data can be found at the publication website as supplementary files.
1) Download the folder with all the scripts to *your-favorite-place* (e.g. Desktop on your local PC). Keep the folder structure as it is. All the files you need are now in *your-favorite-place/DrugInteractionsPipeline/*.
2) Open an RStudio Porject in the folder *your-favorite-place/DrugInteractionsPipeline/*.
3) Download the supplementary files and copy the files 1 to 6 in *your-favorite-place/DrugInteractionsPipeline/InputData*.
4) You are ready to run the pipeline. Start by running *Analysis_of_interactions/Brochado2018_v3.R* and move on to run the scripts inside *Sensitivity Analysis*.

Specific questions can be answered by the authors upon request.

