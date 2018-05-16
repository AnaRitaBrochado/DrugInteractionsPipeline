# DrugInteractionsPipeline
Source code for data analysis - Brochado et al Nature 2018

This folder contains source code implemented and used to analyze the data published by Brochado et al, Nature 2018 (DOI:).
RStudio Version 1.0.136, and R version 2.15.1 were used.
The authors recommend access to a high-performance computing system to run the complete pipeline.
Furthermore, the authors advice caution in re-using this pipeline for other datasets than the one it was built for, due to potential divergences on data quality.

## The pipeline
The pipeline is split in 4 parts:
1) Pre-possessing of raw data & quality control
2) Interaction-scores calculator
3) Analysis of drug interactions
4) Sensitivity analysis

They should be used following the aforementioned order, in order to ensure that all functions and data dependences are satisfied.
The authors recommend a careful look into the directory paths for input-output file allocation in order to ensure that all paths are properly set, before running the entire pipeline.
All figures related to data analysis relevant for the publication mentioned above were generated with the scripts from 3 & 4.

## Input files
The source data can be found at the publication website.
Additional input files needed to run the pipeline are provided here.
The files named “Plate database.txt” are strain dependent and contain all relevant information concerning each multi-well plate, such as a unique plate identifier, batch, run-file, query (also called donor) drug, quality control indicators, etc.
The files named Map.txt are strain dependent and contain all relevant information concerning each well in the multi-well plates, such as array (also called receiver) drug, quality control indicators, etc.

Specific questions can be answered by the authors upon request.

