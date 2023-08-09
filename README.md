# Repository: *Introduced trout alter the trait composition of alpine zooplankton communities*

## Table of Contents
- [Author](#author)
- [Overview](#overview)
- [Contents](#contents)
- [Contents](#contents)
- [Script workflow summary](#script-workflow-summary)
  - [data](#data)
  - [results](#results)
  - [scripts](#scripts)

## Author
Name| Contact
:-----:|:-----
Alex Lewanski|[allewanski AT gmail.com](mailto:allewanski@gmail.com)
Spencer Cruz|[scruz6 AT uwyo.edu](mailto:scruz6@uwyo.edu)
Lindsey Boyle |
Catherine Wagner|
Amy Krist|

## Overview
The `wr_zooplankton_trait_project` repo is associated with the following manuscript: *Introduced trout alter the trait composition of alpine zooplankton communities*. Further details will be added upon acceptance of the manuscript.

## Contents
```
wr_zooplankton_trait_project
|-- data
|   `-- processed_data
|-- scripts
|-- figures
|-- results
|   |-- mod_diagnostics
|   `-- summary_info
```

## Script workflow summary
```
`zooplankton_data_processing_script.R`
| action: initial processing of trait and count data
| inputs: Zooplankton_eggs_and_length_3_29_2022.xlsx ("Data formatted for R" sheet)
|         Rotifer_Measurements_and_Traits.xlsx (trait data for rotifers)
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
|     |        zoop_length_total_list0.rds (processed length data)
|     |        WR_data_ALL_Boyle.xlsx (dataframe of fish presence/absence info for each lake)
|     |        Zooplankton_eggs_and_length_5_9_2022.xlsx ("Inter Traits" sheet)
|     |        all_taxa_trait_info_processed.rds (processed dataframe of species-level trait data for zooplankton taxa)
|     | outputs:
|     |  TABULAR DATA:
|     |     alpha_ses_fd.csv (standardized effect sizes of FRic, FDis, FDiv, and FEve for each lake)
|     |     alpha_scheiner_metrics.csv (Scheiner metric values for each lake)
|     |     fd_summary_table.csv [summaries of all functional diversity metrics (in the `summary_info` directory)]
|     |  PLOTS:
|     |     updated_multidim_FD_violin.png (violin plot of FRic, FDis, FDiv, and FEve in fish vs. fishless lakes)
|     |     scheiner_fd_metrics_boxplot.png (violin plot of Scheiner metrics in fish vs. fishless lakes)
|     |    
|     `--> `alpha_diversity_analyses.R`
|            actions: conducts some initial data exploration (e.g. outliers, collinearity, homogeneity of variance, etc.)
|                     fit models of with each FD variable as a response and the folowing predictors: temp + area + julian_date + time_convert + fish
|                     model diagnostics
|                     model results visualization
|	         notes: script requires custom functions from the following script: `custom_functions_zooplankton_project.R`
|            inputs: WR_data_ALL_Boyle.xlsx (metadata including fish/fishless info, lake chemistry, elevation, etc.)
|             	     alpha_ses_fd.csv (standardized effect sizes of FRic, FDis, FDiv, and FEve)
|           	     alpha_scheiner_metrics.csv (Scheiner metrics)
|            outputs:
|	          TABULAR DATA:
|	             mixed_mod_selection.csv (fit summary for mixed models)
|	          PLOTS:
|	             multidim_fd_alpha.png (boxplots of FRic, FDis, FDiv, and FEve)
|	             scheiner_alpha_plot.png (boxplots of scheiner metrics)
|      
|--> `intraspecific_trait_analyses.R`
|      action: perform intraspecific length analyses
|      notes: script requires custom functions from the following scripts: `custom_functions_zooplankton_project.R` and `intraspecific_variation_functions.R`
|             important sources:
|               [Leps et al. (2011) Ecography](https://doi.org/10.1111/j.1600-0587.2010.06904.x)
|               [de Bello et al. (2011) Methods Eco Evo](https://doi.org/10.1111/j.2041-210X.2010.00071.x)
|      inputs: count_dat_processed_list.rds (processed count data)
|              zoop_length_total_list0.rds (processed length data)
|              WR_data_ALL_Boyle.xlsx (metadata including fish/fishless info, lake chemistry, elevation, etc.)
|              Zooplankton_eggs_and_length_5_9_2022.xlsx ("Inter Traits" sheet)
|              all_taxa_trait_info_processed.rds (processed dataframe of species-level trait data for zooplankton taxa)
|              WR2018ZoopsID.xlsx (all 43 lake sheets)
|              WR2019ZoopsID.xlsx (all 60 lake sheets)             
|      outputs: 
|	    TABULAR DATA:
|         fish_info_anova_processed.csv (results of ANOVA that "decomposes variation of trait values within a community into three sources: (i) the intraspecific trait variability, (ii) the variability due to species turnover and (iii) their covariation is also separated")
|	      vardecomp_fish_vs_fishless.csv [proportion of variance in trait values explained by within vs. between species components (and their covariation)]
|       PLOTS: 
|         intra_turnover_prop_results_treatmentmean_plot.png (barchart of anova showing how much variation in length is explained by fish presence)
|         trait_com_mean_mod_estimate_treatmentmean.png (plot showing community-weighted means in fish vs fishless lakes)
|         species_level_shift_plot.png (plot of length data for each species in fish vs. fishless lakes)
|         vardecomp_fish_vs_fishless_plot.png (barchart of variance explained by w/in vs. between changes in zooplankton length between fish vs. fishless lakes)
|
`--> `zooplankton_beta_fd_analyses.R`
       action: Beta diversity comparisons between fish and fishless lakes
       notes: script requires custom functions from `custom_functions_zooplankton_project.R`
       inputs:
         count_dat_processed_list.rds (processed count data)
         zoop_length_total_list0.rds (processed length data)
         WR_data_ALL_Boyle.xlsx (metadata including fish/fishless info, lake chemistry, elevation, etc.)
         Zooplankton_eggs_and_length_5_9_2022.xlsx ("Inter Traits" sheet)
         all_taxa_trait_info_processed.rds (processed dataframe of species-level trait data for zooplankton taxa)
       outputs: 
 	    TABULAR DATA:
          pairwise_beta_fd_scheiner.csv (beta diversity calculations for all lake pairs)
 	      vardecomp_fish_vs_fishless.csv (PCoA visualization of beta diversity)
        PLOTS: 
          beta_div_pcoa_fig.png [plot of beta diversity values (both weighted means and Scheiner metrics)]   


```


## Content details
### `data`
The main data associated with the project.

### `results`
The results from the project's analyses including metric calculations, outputs of statistical analyses, etc.

### `scripts`
The scripts for analysis, processing, and generation of the figures and tables. All scripts in this directory were used on my local computer.
