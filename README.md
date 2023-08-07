# Repository: *Introduced trout alter the trait composition of alpine zooplankton communities*

## Table of Contents
- [Author](#author)
- [Overview](#overview)
- [Contents](#contents)
- [Content details](#content-details)
  - [data](#data)
  - [results](#results)
  - [scripts](#scripts)

## Author
Name| Contact
:-----:|:-----
Alex Lewanski|[allewanski AT gmail.com](mailto:allewanski@gmail.com)
Spencer Cruz|[email](email@email.com)
Lindsey Boyle | [email](email@email.com)
Catherine Wagner|[email](email@email.com)
Amy Krist|[email](email@email.com)

## Overview
The `wr_zooplankton_trait_project` repo is associated with the following manuscript: *Introduced trout alter the trait composition of alpine zooplankton communities*. Further details will be added upon acceptance of the manuscript.

## Contents
```
wr_zooplankton_trait_project
|-- data
|   `-- processed_data
`-- scripts
```

## Script workflow summary
```
`zooplankton_data_processing_script.R`
| action: initial processing of trait and count data
| inputs: Zooplankton_eggs_and_length_3_29_2022.xlsx ("Data formatted for R" sheet)
|         Rotifer_Measurements_and_Traits.xlsx
|         Zooplankton_eggs_and_length_5_9_2022.xlsx ("Inter Traits" sheet)
|         WR2018ZoopsID.xlsx (all 43 lake sheets)
|         WR2019ZoopsID.xlsx (all 60 lake sheets)
| outputs: all_taxa_trait_info_processed.rds (processed dataframe of species-level trait data for zooplankton taxa)
|          count_dat_processed_list.rds [list of processed count data (split by year)]
|          zoop_length_total_list0.rds [list of processed individual-level length data (split by year)]
|
|--> `zooplankton_FD_analysis.R`
|     | action: calculate the species-level alpha functional diversity calculations
|     | notes: script requires custom functions from the following script: `custom_functions_zooplankton_project.R`
|     | input: count_dat_processed_list.rds (processed count data)
|	  |		 zoop_length_total_list0.rds (processed length data)
|	  |		 WR_data_ALL_Boyle.xlsx (dataframe of fish presence/absence info for each lake)
|     |       Zooplankton_eggs_and_length_5_9_2022.xlsx
|     |       all_taxa_trait_info_processed.rds
|     | outputs:
|	  | TABULAR DATA:
|	  |		alpha_ses_fd.csv (standardized effect sizes of FRic, FDis, FDiv, and FEve for each lake)
|	  |		alpha_scheiner_metrics.csv (Scheiner metric values for each lake)
|	  |		fd_summary_table.csv [summaries of all functional diversity metrics (in the `summary_info` directory)]
|     |	PLOTS:
|	  | 	updated_multidim_FD_violin.png (violin plot of FRic, FDis, FDiv, and FEve in fish vs. fishless lakes)
|	  |	scheiner_fd_metrics_boxplot.png (violin plot of Scheiner metrics in fish vs. fishless lakes)
|     |         
|     |--> `alpha_diversity_analyses.R`
|           actions:
|	  	    notes: script requires custom functions from the following script: `custom_functions_zooplankton_project.R`
|		    inputs: WR_data_ALL_Boyle.xlsx (metadata including fish/fishless info, lake chemistry, elevation, etc.)
|				   alpha_ses_fd.csv (standardized effect sizes of FRic, FDis, FDiv, and FEve)
|				   alpha_scheiner_metrics.csv (Scheiner metrics)
|



```


## Content details
### `data`
The main data associated with the project.

### `results`
The results from the project's analyses including metric calculations, outputs of statistical analyses, etc.

### `scripts`
The scripts for analysis, processing, and generation of the figures and tables. All scripts in this directory were used on my local computer.
