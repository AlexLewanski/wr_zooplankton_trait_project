#################################################
#################################################
### ANALYSES OF INTRASPECIFIC TRAIT VARIATION ###
#################################################
#################################################


##########################
### SCRIPT PREPARATION ###
##########################

### PACKAGES ###
library(here)
library(tidyverse)
library(gawdis)
library(readxl)
library(picante) #matrix randomization
library(ggtext)
library(ggridges)

#load custom functions
source(here("scripts", "custom_functions_zooplankton_project.R"))
source(here("scripts", "intraspecific_variation_functions.R"))


### LOADING DATA ###
#(processed) count and length data
count_dat_processed_list <- readRDS(file = here('data', 'processed_data', 'count_dat_processed_list.rds'))
zoop_length_total_list0 <- readRDS(file = here('data', 'processed_data', 'zoop_length_total_list0.rds'))

#fish presence/absence info
lindsey_master_data <- as.data.frame(read_excel(here('data', 'WR_data_ALL_Boyle.xlsx')))

#trait data (not including length data)
trait_info_init <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_5_9_2022.xlsx'), sheet = 'Inter Traits'))

all_taxa_traits <- readRDS(here('data', 'processed_data', 'all_taxa_trait_info_processed.rds'))



##########################
### PROCESSING OF DATA ###
##########################

#*** most processing done in zooplankton_data_processing_script.R

#focal taxa included in analyses
focal_taxa <- c('b. longirostris',
                'd. mendatoe',
                'd. middenorfiana/pulex',
                "h. shoshone",
                'l. minutus',
                'chydorus',
                'h. gibberum',
                'p. pediculus')


###  FISH PRESENCE/ABSENCE INFO FOR EACH LAKE ###
fish_info <- lindsey_master_data %>% 
  select(lake, fish)

#fix lake names in the fish presence/absence dataframe
processed_lake_info <- fish_info[!duplicated(fish_info$lake),] %>%
  mutate(lake = recode(lake,
                       'lake37' = "lake_37",
                       'blackjoe' = "black_joe",
                       'upper_blackjoe' = "upper_black_joe",
                       "below_tb" = "lake_below_tb"))


###  TRAIT DATA ###
#remove euchlanis and filinia and move species to a row name
all_taxa_traits_processed_init <- all_taxa_traits %>% 
  filter(!(species %in% c('euchlanis', 'filinia')) ) %>%
  column_to_rownames(var = 'species')


###  ABUNDANCE DATA ###

count_dat_processed_init <- count_dat_processed_list %>%
  bind_rows() %>% 
  filter(Lake != "Lost Lake") %>% #remove lost lake
  mutate(Lake = dplyr::recode(Lake, "N of Blue" = "North Blue")) %>% #fix lake name
  group_by(Lake) %>% 
  filter(year == max(year)) %>% #keep only the most recent survey
  ungroup() %>% 
  group_by(Lake, taxa) %>%
  summarize(total_abund = sum(total), .groups = 'drop') %>% 
  filter(taxa %in% unique(rownames(all_taxa_traits_processed_init) )) %>%
  mutate(lake = gsub(" ", "_", str_trim(tolower(Lake)))) %>%
  mutate(lake = dplyr::recode(lake, 
                              'annika_lake' = 'annika',
                              "bear_west" = "west_bear",
                              "east_no_name" = "east_noname",
                              "frozen_lake_1" = "frozen1",
                              "frozen_lake_2" = "frozen2",
                              "island_lake" = 'island',
                              "lindsey_lake" = "lindsey",
                              "little_mountain_sheep" = "little_mtnsheep",
                              "low_deep_creek" = "lower_deepcreek",
                              "lower_indian_basin" = "lower_indian",
                              "mid_deep_creek" = "mid_deepcreek",
                              "mid_indian_basin" = "mid_indian",
                              "n_of_blue" = "north_blue",
                              "no_name_west" = "west_noname",
                              "upper_indian_basin" = "upper_indian",
                              "west_of_mt_lester" = "west_mtlester")) %>% 
  select(!Lake)

all_taxa_traits_processed <- all_taxa_traits_processed_init %>% 
  filter(rownames(all_taxa_traits_processed_init) %in% unique(count_dat_processed_init$taxa))


zoop_length_total_list0 <- list()

for (YEAR in c(2018, 2019)) {
  
  wr_zoop_init <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
    message(NAME)
    len_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 5, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
    
    #remove all extra columns (everything that is found to the right of the first column labelled as ...NUMBER)
    #len_dat_1 <- len_dat[,!grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))]
    extra_col_detect <- grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))
    if (any(extra_col_detect)) {
      remove_col_indices <- min(which(extra_col_detect)):ncol(len_dat)
      len_dat_1 <- len_dat[, !(seq_len(ncol(len_dat)) %in% remove_col_indices) ]
    } else {
      len_dat_1 <- len_dat
    }
    
    len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(grepl(pattern = 'Averages', len_dat_1[,1])):nrow(len_dat_1)),] #remove the average information
    #len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(len_dat_1[,1] == 'Averages'):nrow(len_dat_1)),]
    
    len_dat_3 <- len_dat_2[apply(len_dat_2, MARGIN = 1, function(x) !all(is.na(x))),] #remove any rows with all NAs (e.g. extra rows at the bottom of the dataframe)
    colnames(len_dat_3) <- tolower(colnames(len_dat_3)) #sometimes "1x" is labelled as "1X", which messes things up
    
    return(len_dat_3)
  })
  
  
  wr_zoop_1 <- lapply(wr_zoop_init, function(DATA) {
    as.data.frame(apply(DATA, 2, as.numeric)) %>%
      mutate(id = 1:n()) %>% #add unique id for each row (necessary for correct wide pivot)
      pivot_longer(cols = !id, names_to = 'taxa_trait', values_to = 'value', values_drop_na = TRUE) %>% 
      mutate(taxa = tolower(str_trim(str_remove(string = taxa_trait, pattern = "1x$|2x$|mm$"))),
             taxa = ifelse(taxa %in% c("polyphemus"), "p. pediculus", taxa),
             taxa = gsub("\\.$", "", taxa),
             measurement = paste0('length_', tolower(str_trim(str_extract(string = taxa_trait, pattern = "1x$|2x$|mm$"))))) %>% 
      select(!taxa_trait) %>% 
      pivot_wider(names_from = measurement, values_from = value) %>% 
      select(!id)
  })
  
  zoop_length_total_list0[[paste0('length_dat_', YEAR)]] <- as.data.frame(bind_rows(wr_zoop_1, .id = 'Lake'))
  zoop_length_total_list0[[paste0('length_dat_', YEAR)]]$year <- YEAR
  zoop_length_total_list0[[paste0('length_dat_', YEAR)]] <- zoop_length_total_list0[[paste0('length_dat_', YEAR)]] %>% 
    mutate(Lake = recode(Lake, "north_blue" = "n of blue"))
  
}


survey_year <- bind_rows(zoop_length_total_list0) %>% 
  mutate(year_lake = paste0(Lake, year)) %>% 
  group_by(Lake) %>% 
  filter(year == max(year)) %>% 
  slice_head(n = 1) %>% 
  pull(year_lake)

subset_length_year <- zoop_length_total_list0 %>% 
  bind_rows() %>% 
  mutate(year_lake = paste0(Lake, year)) %>% 
  filter(year_lake %in% survey_year)

length_processed <- subset_length_year %>%
  mutate(taxa = recode(taxa, 
                       'B_longirostris' = 'b. longirostris',
                       'D_mendotae' = 'd. mendatoe',
                       'D_pulex' = 'd. middenorfiana/pulex',
                       'H_shoshone' = "h. shoshone",
                       'L_minutus' = 'l. minutus',
                       'Chydorus' = 'chydorus',
                       'H. gibberum' = 'h. gibberum',
                       'P. pediculous' = 'p. pediculus')) %>% 
  filter(taxa != 'harpacticoida') %>%
  filter(Lake != "Lost Lake") %>% #remove lost lake
  mutate(lake = gsub(" ", "_", str_trim(tolower(Lake)))) %>%
  mutate(lake = recode(lake,
                       'annika_lake' = 'annika',
                       "bear_west" = "west_bear",
                       "east_no_name" = "east_noname",
                       "frozen_lake_1" = "frozen1",
                       "frozen_lake_2" = "frozen2",
                       "island_lake" = 'island',
                       "lindsey_lake" = "lindsey",
                       "little_mountain_sheep" = "little_mtnsheep",
                       "low_deep_creek" = "lower_deepcreek",
                       "lower_indian_basin" = "lower_indian",
                       "mid_deep_creek" = "mid_deepcreek",
                       "mid_indian_basin" = "mid_indian",
                       "n_of_blue" = "north_blue",
                       "no_name_west" = "west_noname",
                       "upper_indian_basin" = "upper_indian",
                       "west_of_mt_lester" = "west_mtlester",
                       "bap_lake" = "bap",
                       "cutthroat_(stough)" = "cutthroat_stough")) %>% 
  filter(taxa %in% focal_taxa) %>% 
  mutate(taxa_lake = paste0(taxa, "SPACE", lake))

count_dat_processed1 <- count_dat_processed_init %>%
  mutate(lake = recode(lake,
                       "bap_lake" = "bap",
                       "cutthroat_(stough)" = "cutthroat_stough")) %>% 
  filter(taxa %in% focal_taxa) %>% 
  group_by(lake) %>% 
  mutate(rel_abund = total_abund/sum(total_abund),
         taxa_lake = paste0(taxa, "SPACE", lake)) %>% 
  ungroup()


### IDENTIFYING PROBLEM LAKES (CONTAIN TAXA WITH COUNT DATA BUT NOT LENGTH DATA)###
#*** THIS SHOULD BE BETTER RESOLVED IN THE FINAL ANALYSIS
taxa_count_not_length <- unique(count_dat_processed1$taxa_lake)[!unique(count_dat_processed1$taxa_lake) %in% length_processed$taxa_lake]
unique(length_processed$taxa_lake)[!unique(length_processed$taxa_lake) %in% unique(count_dat_processed1$taxa_lake)]

#lakes with count data but no length data (remove these)
problem_lakes <- unique(gsub('.*SPACE', '', taxa_count_not_length))
count_dat_processed1_subset <- count_dat_processed1 %>% 
  filter(!lake %in% problem_lakes)



################################################################################
### INTRASPECIFIC SHIFTS IN COMPOSITION (TREATMENT SUMMARIES OF MEAN LENGTH) ###
################################################################################

#community means for each taxa in the fish and fishless lakes
length_mean_info_fish_init <- left_join(length_processed, processed_lake_info, by = 'lake') %>% 
  mutate(fish = as.character(fish))

#size of each taxa in the fish vs. fishless lakes
length_mean_info_fish <- length_mean_info_fish_init %>% 
  group_by(taxa) %>% 
  mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>%
  ungroup() %>% 
  group_by(fish, taxa) %>% 
  summarize(sample_size = n(),
            fish_mean_length = mean(length_mm, na.rm = TRUE),
            overall_mean_length = first(overall_mean_length),
            .groups = 'drop') %>% 
  mutate(taxa_fish = paste(taxa, fish, sep = "_"))

# community level (i.e., lake level) means of zooplankton length
processed_mean_length_fish <- left_join(count_dat_processed1_subset, processed_lake_info, by = 'lake') %>% 
  mutate(taxa_fish = paste(taxa, fish, sep = "_")) %>%
  left_join(., length_mean_info_fish, by = 'taxa_fish') %>% 
  group_by(lake) %>% 
  summarize(com_specific_length = weighted.mean(fish_mean_length, w = rel_abund),
            overall_length = weighted.mean(overall_mean_length, w = rel_abund),
            fish = as.factor(first(fish.x)),
            .groups = 'drop')

#trait flex anova based on zooplankton community length
fish_info_anova <- trait.flex.anova(~fish, com_specific_length, overall_length, 
                                    data = processed_mean_length_fish)

fish_info_anova_processed <- fish_info_anova$RelSumSq %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'component') %>%
  mutate(component = trimws(component)) %>%
  pivot_longer(cols = !component, names_to = 'component2', values_to = 'proportion_variance')

intra_turnover_prop_results_treatmentmean_plot <- fish_info_anova_processed
  mutate(percentage_variance = proportion_variance*100,
         trait = 'length',
         component_processed = factor(case_when(component == 'fish' ~ 'Fish presence',
                                                component == 'Residuals' ~ 'Error',
                                                component == 'Total' ~ 'Total'),
                                      levels = c('Fish presence', 'Error', 'Total')),
         component2_processed = factor(case_when(component2 == 'Covariation' ~ 'Covariation',
                                                 component2 == 'Intraspec.' ~ 'Intraspecific shift',
                                                 component2 == 'Turnover' ~ 'Interspecific shift'),
                                       levels = c('Covariation', 'Intraspecific shift', 'Interspecific shift'))) %>% 
  ggplot() +
  geom_bar(data = . %>% filter(component2 != 'Total'),
           mapping = aes(x = trait, y = percentage_variance, fill = component2_processed), 
           stat = 'identity', position = 'stack') +
  facet_wrap(~component_processed) +
  ylab('% explained variance') +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "#d8d8d8", linewidth = 1, color = "#d8d8d8"),
        panel.border = element_rect(fill = NA, linewidth = 1, color = "#d8d8d8"),
        strip.text = element_text(size = 14, face = "plain", color = "black")) +
  theme(axis.title.x = element_blank())

write.csv(fish_info_anova_processed, 
          here('results', 'summary_info', 'fish_info_anova_processed.csv'),
          quote = FALSE,
          row.names = FALSE)
  
ggsave(plot = intra_turnover_prop_results_treatmentmean_plot,
       filename = here('figures', 
                       'intraspecific_results', 
                       'intra_turnover_prop_results_treatmentmean_plot.png'), 
       width = 15*0.7, height = 8*0.7, device = 'png')

#ggsave(plot = intra_turnover_prop_results_treatmentmean_plot,
#       filename = here('figures', 'intraspecific_results', 'intra_turnover_prop_results_treatmentmean_plot.pdf'), 
#       width = 15*0.7, height = 8*0.7, device = 'pdf')

confidence_intervals_fish <- lapply(setNames(nm = c('com_specific_length', 'overall_length')), function(x, dat) {
  predict(lm(as.formula(paste0(x, ' ~ fish')), data = dat), 
          data.frame(fish = factor(c('0', '1')) ),
          interval="confidence", level = 0.95) %>% 
    as.data.frame() %>% 
    mutate(outcome_var = x,
           fish_info = c('Fishless', 'Fish'))
}, dat = processed_mean_length_fish) %>% 
  bind_rows() %>% 
  `rownames<-`( NULL )


trait_com_mean_mod_estimate_treatmentmean_plot <- confidence_intervals_fish %>%
  mutate(outcome_var_final = if_else(outcome_var == 'com_specific_length', 'Specific', 'Fixed')) %>% 
  ggplot() +
  geom_point(data = processed_mean_length_fish %>% 
                pivot_longer(cols = c(com_specific_length, overall_length), names_to = 'outcome_var', values_to = 'fit') %>% 
                mutate(fish_info = if_else(fish == '0', 'Fishless', 'Fish')) %>%
               mutate(outcome_var_final = if_else(outcome_var == 'com_specific_length', 'Specific', 'Fixed')),
              aes(x = fish_info, y = fit, color = outcome_var_final),
              position = position_jitterdodge(dodge.width = 0.75,  jitter.width = 0.1),
              size = 2.25, alpha = 0.3) +
  geom_line(aes(x = fish_info, y = fit, group = outcome_var_final, color = outcome_var_final),
            position = position_dodge(width = 0.75), linewidth = 2) +
  geom_errorbar(aes(x = fish_info, ymin = lwr, ymax = upr, color = outcome_var_final),
                position = position_dodge(width = 0.75), width = 0.08, linewidth = 2) +
  geom_point(aes(x = fish_info, y = fit, color = outcome_var_final),
             position = position_dodge(width = 0.75),
             size = 7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab('Length')

ggsave(plot = trait_com_mean_mod_estimate_treatmentmean_plot,
       filename = here('figures', 'intraspecific_results', 'trait_com_mean_mod_estimate_treatmentmean.png'), 
       width = 11*0.7, height = 8*0.7, device = 'png')

#ggsave(plot = trait_com_mean_mod_estimate_treatmentmean_plot,
#       filename = here('figures', 'intraspecific_results', 'trait_com_mean_mod_estimate_treatmentmean.pdf'), 
#       width = 11*0.7, height = 8*0.7, device = 'pdf')

#filter dataset to only the taxa found in both fish and fishless lakes
taxa_inboth_treatments <- length_mean_info_fish_init %>%
  filter(!is.na(length_mm)) %>% 
  group_by(taxa) %>% 
  filter(all(c('1', '0') %in% fish)) %>% 
  ungroup() %>% 
  pull(taxa) %>% 
  unique()

#order taxa names in descending order (as a vector)
taxa_descending_order <- length_mean_info_fish_init %>% 
  filter(!is.na(length_mm)) %>%
  group_by(taxa) %>% 
  #filter(all(c('1', '0') %in% fish)) %>%
  summarize(mean_val = mean(length_mm)) %>% 
  arrange(desc(mean_val)) %>% 
  pull(taxa)

#conduct Mann–Whitney U test on length data for fish vs. fishless lakes
mannu_results <- lapply(setNames(nm = taxa_inboth_treatments), function(x, dat) {
  mannu_output <- wilcox.test(length_mm ~ fish, 
              data = dat %>% filter(taxa == x), 
              na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
  
  return(data.frame(taxa = x,
                    W = mannu_output$statistic,
                    P = ifelse(mannu_output$p.value < 0.001, '< 0.001', paste0('= ', round(mannu_output$p.value, 3)) ))
         )
  
}, dat = length_mean_info_fish_init) %>% 
  bind_rows() %>% 
  mutate(factor = factor(taxa, levels = taxa_descending_order),
         results_text = paste0('W = ', W, '<br>P ', P))

#calculate mean length of fishless and fish lakes and the overall mean
length_info_processed <- length_mean_info_fish %>% 
  #group_by(taxa) %>% 
  #filter(all(c('1', '0') %in% fish)) %>%
  #ungroup() %>% 
  mutate(fish_info = factor(if_else(fish == '0', 'Fishless', 'Fish'), levels = c('Fishless', 'Fish')),
         taxa = factor(taxa, levels = taxa_descending_order)) %>% 
  select(fish_info, taxa, fish_mean_length, overall_mean_length) %>% 
  pivot_wider(names_from = 'fish_info', values_from = 'fish_mean_length') %>% 
  dplyr::rename(`Mean (fishless)` = Fishless,
                `Mean (fish)` = Fish,
                `Fixed mean` = overall_mean_length) %>% 
  pivot_longer(cols = !taxa, names_to = 'mean_type', values_to = 'value') %>% 
  as.data.frame()

#plot of intraspecific lengths
species_level_shift_plot <- length_mean_info_fish_init %>%
  filter(!is.na(length_mm)) %>% 
  #group_by(taxa) %>% 
  #filter(all(c('1', '0') %in% fish)) %>%
  #ungroup() %>% 
  mutate(fish_info = factor(if_else(fish == '0', 'Fishless', 'Fish'), levels = c('Fishless', 'Fish')),
         taxa = factor(taxa, levels = taxa_descending_order)) %>% 
  ggplot(aes(y = taxa, x = length_mm, fill = fish_info, color = fish_info)) + #height = ..count..
  geom_density_ridges(stat = "density_ridges",
                      scale = 0.98, rel_min_height = .01, alpha = 0.5,
                      point_size = 0.5, size = 0.25, #point_shape = "|",
                      position = position_points_jitter(width = 0.03, height = 0.05, yoffset = -0.03) ) +
  geom_jitter(data = . %>% 
                filter(fish_info == 'Fishless'),
              #aes(y = as.numeric(taxa) - 0.03),
              aes(y = as.numeric(taxa) - 0.03),
              size = 0.3, alpha = 0.7,
              position = position_jitterdodge(dodge.width = 0,
                                              jitter.width = 0.07,
                                              jitter.height = 0)) +
  geom_jitter(data = . %>% filter(fish_info == 'Fish'),
              aes(y = as.numeric(taxa) - 0.08),
              size = 0.3, alpha = 0.7,
              position = position_jitterdodge(dodge.width = 0,
                                              jitter.width = 0.07,
                                              jitter.height = 0)) +
  geom_point(data = length_info_processed %>% filter(mean_type == 'Fixed mean'),
             mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#495057', fill = '#495057') +
  geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fishless)' & !is.na(value)),
             mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#c1121f', fill = '#c1121f') +
  geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fish)' & !is.na(value)),
             mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#662e9b', fill = '#662e9b') +
  geom_textbox(
    data = mannu_results,
    aes(y = taxa, x = Inf, label = results_text),
    hjust = 1, halign = 0, size = 3,
    box.colour = NA, fill = NA, # Hide the box fill and outline
    box.padding = unit(rep(2.75, 4), "pt"), colour = "black",
    vjust = 1, nudge_y = 0.37, width = NULL
  )  +
  #annotate("text", x = 0, y = as.numeric(prac_df$taxa) + 0.5, label = prac_df$test) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  xlab('Length (mm)')


ggsave(plot = species_level_shift_plot,
       filename = here('figures', 'intraspecific_results', 'species_level_shift_plot.png'), 
       width = 10*0.7, height = 11*0.7, device = 'png')



#################################################################################
### INTRASPECIFIC SHIFTS IN COMPOSITION (LAKE-LEVEL SUMMARIES OF MEAN LENGTH) ###
#################################################################################

count_length_processed1_lake <- left_join(count_dat_processed1_subset, 
                                     length_processed %>% 
                                       group_by(taxa) %>% 
                                       mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>% 
                                       ungroup() %>% 
                                       group_by(taxa_lake) %>% 
                                       summarize(sample_size = n(),
                                                 com_mean_length = mean(length_mm, na.rm = TRUE),
                                                 overall_mean_length = first(overall_mean_length),
                                                 .groups = 'drop'), 
                                     by = 'taxa_lake') %>%
  filter(!is.nan(overall_mean_length) | !is.nan(com_mean_length)) 

count_length_processed2_lake <- count_length_processed1_lake %>%
  group_by(lake) %>% 
  summarize(com_specific_length = weighted.mean(com_mean_length, w = rel_abund),
            overall_length = weighted.mean(overall_mean_length, w = rel_abund), 
            .groups = 'drop')

count_length_processed2 <- left_join(count_length_processed2_lake,
                                     processed_lake_info,
                                     by = 'lake') %>% 
  mutate(fish = factor(fish))


lake_level_anova <- trait.flex.anova(~fish, com_specific_length, overall_length, 
                                     data = count_length_processed2)


lake_level_anova$RelSumSq %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'component') %>%
  mutate(component = trimws(component)) %>%
  pivot_longer(cols = !component, names_to = 'component2', values_to = 'proportion_variance') %>% 
  #filter(component2 != 'Covariation') %>%
  mutate(percentage_variance = proportion_variance*100,
         trait = 'length',
         component_processed = factor(case_when(component == 'fish' ~ 'Fish presence',
                                                component == 'Residuals' ~ 'Error',
                                                component == 'Total' ~ 'Total'),
                                      levels = c('Fish presence', 'Error', 'Total')) ) %>% 
  ggplot() +
  geom_bar(data = . %>% filter(component2 != 'Total'),
           mapping = aes(x = trait, y = percentage_variance, fill = component2), 
           stat = 'identity', position = 'stack') +
  facet_wrap(~component_processed) +
  ylab('% explained variance') +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


#exploration of sample sizes underlying mean values
count_length_processed1_lake$sample_size

count_length_processed1_lake %>%
  ggplot() +
  geom_histogram(aes(x = sample_size), bins = 30) +
  theme_bw()



##############################
### VARIANCE DECOMPOSITION ###
##############################

#goal: calculate proportion of variance in trait values explained by within vs. between
#species components (and their covariation)

#important paper: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2010.00071.x

#As a first attempt, I am quantifying the proportion within the fish vs. fishless lakes and
#I am using the average relative abundances for each species across lakes in these calculations

#step 1: Add zeros in for each species in the lakes that they are absent
count_dat_processed1_addzeros <- count_dat_processed1 %>% 
  select(taxa, lake, rel_abund) %>% 
  pivot_wider(values_from = rel_abund, names_from = 'taxa', values_fill = 0) %>% 
  pivot_longer(!lake, names_to = 'taxa', values_to = 'rel_abund')


#step 2: merge count data with fish info and then calculate average abundance within the fish
#and fishless lakes
combined_rel_abund_fishinfo <- left_join(count_dat_processed1_addzeros,
          processed_lake_info,
          by = 'lake') %>% 
  mutate(fish = as.character(fish)) %>% 
  group_by(fish, taxa) %>% 
  summarize(rel_abund = mean(rel_abund), .groups = 'drop') %>% 
  filter(rel_abund > 0) %>% 
  rename(species = taxa)

#split relative abundance dataframe into list of dataframes based on fish presence
rel_abund_split <- combined_rel_abund_fishinfo %>% 
  group_by(fish) %>% 
  group_split(.keep = FALSE)
names(rel_abund_split) <- unique(combined_rel_abund_fishinfo$fish)

#step 3: merge length data with fish info and split into fish and fishles dataframes
length_fish_split <- left_join(length_processed, processed_lake_info, by = 'lake') %>% 
  select(fish, taxa, length_mm) %>%
  rename(species = taxa, trait = length_mm) %>% 
  filter(!is.na(trait)) %>% 
  group_by(fish) %>% 
  group_split(.keep = FALSE)
names(length_fish_split) <- unique(left_join(length_processed, processed_lake_info, by = 'lake')$fish)

#step 4: create 1000 bootstrap datasets for the length data by resampling (with replacement) the length 
#data for each species within the fish and fishless lake datasets
bootstrap_length_list <- lapply(length_fish_split, function(x) {
  replicate(1000, x %>% 
              group_by(species) %>% 
              sample_n(size = n(), replace = TRUE) %>%
              ungroup(), simplify = FALSE)
})

#step 5: decompose variance within the fish and fishless datasets based on the length data and
#relative abundance information
observed_vardecomp <- lapply(setNames(nm = c('0', '1')), function(x, rel_abund_list, length_list) {
  var_decomp(abundance_df = rel_abund_list[[x]], 
             trait_df = length_list[[x]]) %>% 
    mutate(fish_info = switch(x, "0" = 'fishless', "1" = 'fish'))
  
}, rel_abund_list = rel_abund_split, length_list = length_fish_split) %>% 
  bind_rows()



#step 6: estimate confidence intervals on the variance decomposition values using the bootstrapping
bs_dataset_list <- list()
ci_vardecomp_list <- list()
for (i in names(bootstrap_length_list)) {
  fish_label <- switch(i, "0" = 'fishless', "1" = 'fish')
  
  bs_dataset_list[[fish_label]] <- lapply(bootstrap_length_list[[i]], function(x, rel_abund_df) {
    var_decomp(abundance_df = rel_abund_df, 
               trait_df = x)
    
  }, rel_abund_df = rel_abund_split[[i]]) %>% 
    bind_rows()
  
  ci_vardecomp_list[[fish_label]] <- bs_dataset_list[[fish_label]] %>% 
    pivot_longer(cols = everything(), names_to = 'quantity', values_to = 'value') %>% 
    group_by(quantity) %>% 
    summarize(lower_ci = quantile(value, probs = c(0.025)),
              upper_ci = quantile(value, probs = c(0.975)))
  
  message('finished ', fish_label)
}


#step 7: process variance decomp results
ci_vardecomp_df <- ci_vardecomp_list %>% 
  bind_rows(.id = 'fish_info') #%>% 
  #pivot_wider(names_from = 'quantity', values_from = ends_with('ci') )

observed_vardecomp_pivot <- pivot_longer(observed_vardecomp, cols = !fish_info, names_to = 'quantity', values_to = 'observed')
vardecomp_fish_vs_fishless <- full_join(observed_vardecomp_pivot, ci_vardecomp_df, by = c('fish_info', 'quantity'))

vardecomp_fish_vs_fishless_plot <- vardecomp_fish_vs_fishless %>%
  mutate(comparison_type = gsub("_.*", "", quantity)) %>% 
  filter(quantity %in% c('within_var', 'between_var')) %>% 
  mutate(quantity1 = factor(if_else(quantity == 'within_var', 'within species', 'between species'), levels = c('within species', 'between species')),
         fish_info_update = factor(if_else(fish_info == 'fishless', "Fishless", "Fish"), levels = c('Fishless', 'Fish'))) %>% 
  group_by(fish_info) %>% 
  mutate(prop = observed/sum(observed),
         prop_text = stringi::stri_pad_right(round(prop, 2), 4, 0)) %>% 
  ggplot(mapping = aes(x = fish_info_update, y = observed, fill = quantity1)) +
  geom_bar(data = . %>% filter(quantity %in% c('within_var', 'between_var')),
           stat = 'identity', position = position_dodge()) + #stack
  geom_errorbar(data = . %>% filter(quantity %in% c('within_var', 'between_var')),
                mapping = aes(x = fish_info_update, ymin = lower_ci, ymax = upper_ci),
                position = position_dodge(0.9),
                width = 0.1) +
  geom_point(data = . %>% filter(quantity %in% c('within_var', 'between_var')),
             mapping = aes(x = fish_info_update, y = observed), 
             position = position_dodge(0.9),
             size = 2) +
  ylab('Variance')  +
  geom_text(aes(y = 0.025, label = prop_text), position = position_dodge(0.9), 
            vjust = 0, size = 5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

write.csv(vardecomp_fish_vs_fishless, 
          here('results', 'summary_info', 'vardecomp_fish_vs_fishless.csv'),
          quote = FALSE,
          row.names = FALSE)

ggsave(plot = vardecomp_fish_vs_fishless_plot,
       filename = here('figures', 'intraspecific_results', 'vardecomp_fish_vs_fishless_plot.png'), 
       width = 15*0.7, height = 8*0.7, device = 'png')



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

#visualize results
# vardecomp_fish_vs_fishless_plot <- vardecomp_fish_vs_fishless %>%
#   filter(quantity %in% c('within_var_prop', 'between_var_prop')) %>% 
#   mutate(quantity1 = factor(if_else(quantity == 'within_var_prop', 'within species', 'between species'), levels = c('within species', 'between species'))) %>% 
#   ggplot() +
#   geom_bar(aes(x = fish_info, y = observed, fill = quantity1), stat = 'identity', position = 'stack') +
#   geom_errorbar(data = . %>% filter(quantity1 == 'between species'), #geom_linerange
#     mapping = aes(x = fish_info, ymin = lower_ci, ymax = upper_ci),
#                 width = 0.1) +
#   geom_point(data = . %>% filter(quantity1 == 'between species'),
#                 mapping = aes(x = fish_info, y = observed), size = 4) +
#   ylim(0, 1) +
#   ylab('Variance') +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.title.x = element_blank())
# 
# ggsave(plot = vardecomp_fish_vs_fishless_plot,
#        filename = here('figures', 'intraspecific_results', 'vardecomp_fish_vs_fishless_plot.png'), 
#        width = 15*0.7, height = 8*0.7, device = 'png')


# confint(aov(com_specific_length ~ fish, data = processed_mean_length_fish),
#         level = 0.95) %>% 
#   as.data.frame() %>% 
#   rename(lower_ci = `2.5 %`,
#          upper_ci = `97.5 %`) %>% 
#   mutate(fish_info = factor(c('fish', 'fishless'), levels = c('fishless', 'fish'))) %>% 
#   ggplot() +
#   geom_errorbar(mapping = aes(x = fish_info, ymin = lower_ci, ymax = upper_ci),
#                 width = 0.1) +
#   geom_point(data = estimated_vals, aes(x = fish_info, y = estimate) )
# 
# estimated_vals <- data.frame(fish_info = c('fish', 'fishless'),
#                              estimate = c(coef(lm(com_specific_length ~ fish, data = processed_mean_length_fish))[1], coef(lm(com_specific_length ~ fish, data = processed_mean_length_fish))[2])
#                              

# tind1<-c(1:6)
# tind2<-c(2:9)
# tind3<-c(1:10)
# tind4<-c(c(tind2*2), c(tind2*2)-2)
# 
# trait_info <- data.frame(species = c(rep('sp1', 6), rep('sp2', 8), rep('sp3', 10), rep('sp4', 16)),
#                          trait = c(c(1:6), c(2:9), c(1:10), c(c(tind2*2), c(tind2*2)-2)) 
# )
# 
# 
# trait_info_alt <- data.frame(species = c(rep('sp1', 1), rep('sp2', 8), rep('sp3', 10), rep('sp4', 16)),
#                              trait = c(c(1), c(2:9), c(1:10), c(c(tind2*2), c(tind2*2)-2)) 
# )
# 
# 
# # abund_info <- stack(table(trait_info$species)) %>% 
# #   select(ind, values) %>% 
# #   mutate(ind = as.character(ind)) %>%
# #   rename(species = ind,
# #          abund = values)
# 
# trait_info1 <- trait_info %>% 
#   group_by(species) %>% 
#   mutate(rel_freq = 1/n()) %>% 
#   ungroup()
# 
# 
# abund_info <- data.frame(species = paste0('sp', 1:4),
#                          rel_abund = c(0.1, 0.2, 0.3, 0.4))
# 
# var_decomp(abundance_df = abund_info, trait_df = trait_info_alt[sample(nrow(trait_info_alt)),])
# 
# var_decomp(abundance_df = trait_info %>% 
#              group_by(species) %>% 
#              summarize(abund = n()) %>% 
#              ungroup() %>% 
#              mutate(rel_abund = abund/nrow(trait_info)) %>% 
#              select(species, rel_abund), 
#            trait_df = trait_info)
# 
# 
# length_processed
# 
#
#
# taxa_count_not_length <- unique(count_dat_processed1$taxa_lake)[!unique(count_dat_processed1$taxa_lake) %in% length_processed$taxa_lake]
# unique(length_processed$taxa_lake)[!unique(length_processed$taxa_lake) %in% unique(count_dat_processed1$taxa_lake)]
# 
# problem_lakes <- unique(gsub('.*SPACE', '', taxa_count_not_length))
# 
# 
# 
# length_processed %>% 
#   filter(lake == 'cliff') %>% 
#   group_by(taxa) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   select(lake, taxa, taxa_lake)
# 
# count_dat_processed1 %>% 
#   filter(lake == 'cliff') %>% 
#   select(lake, taxa, taxa_lake)
# 
# 
# unique(length_processed$taxa_lake)[unique(length_processed) %in% count_dat_processed1$taxa]
# 
# length_processed %>% 
#   filter(lake == '47aaron') %>% 
#   group_by(taxa) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   select(Lake, taxa, taxa_lake)
#   
# 
# count_dat_processed1_subset <- count_dat_processed1 %>% 
#   filter(!lake %in% problem_lakes)
# 
# 
# count_length_processed1 <- left_join(count_dat_processed1_subset, 
#           length_processed %>% 
#             group_by(taxa) %>% 
#             mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>% 
#            ungroup() %>% 
#             group_by(taxa_lake) %>% 
#             summarize(com_mean_length = mean(length_mm, na.rm = TRUE),
#                       overall_mean_length = first(overall_mean_length),
#                       .groups = 'drop'), 
#           by = 'taxa_lake') %>%
#   filter(!is.nan(overall_mean_length) | !is.nan(com_mean_length)) %>%
#   group_by(lake) %>% 
#   summarize(com_specific_length = weighted.mean(com_mean_length, w = rel_abund),
#             overall_length = weighted.mean(overall_mean_length, w = rel_abund), 
#             .groups = 'drop')
# 
# count_length_processed2 <- left_join(count_length_processed1,
#                                      processed_lake_info,
#                                      by = 'lake') %>% 
#   mutate(fish = factor(fish))
# 
# unique(length_processed$taxa_lake)[!unique(length_processed$taxa_lake) %in% count_dat_processed1$taxa_lake]
# unique(count_dat_processed1$taxa_lake)[unique(count_dat_processed1$taxa_lake) %in% length_processed$taxa_lake]
# 
# 
# unique(length_processed$lake)[!(unique(length_processed$lake) %in% count_dat_processed_init$lake)]
# 
# unique(length_processed$taxa)[!unique(length_processed$taxa) %in% unique(count_dat_processed1$taxa)]
# 
# 
# 
# count_dat_processed_init %>% 
#   group_by(lake) %>% 
#   mutate(rel_abund = total_abund/sum(total_abund))
# 
# length_processed %>% 
#   mutate
# 
# unique(length_processed$Lake[!(length_processed$Lake %in% processed_lake_info$lake)])
# 
# 
# 
# left_join(length_processed, processed_lake_info, by = 'lake') %>%
#   group_by(taxa) %>% 
#   mutate(species_mean = mean(length_mm))
#   mutate(fish = factor(fish)) %>% 
#   select(fish, length_mm)
#   
#   
# table(length_processed$lake,
#       length_processed$taxa)
# 
# 
# 
# source('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/STOICH/wind_river_zooplankton/scripts/e6904_trait-flex-v3.r')
# 
# 
# height.trait<- data.frame(MOWING=as.factor(c(0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0))) 
# height.trait$FERTIL <- as.factor(c(0,1,0,1,1,0,1,0,1,0,1,0))
# height.trait$specific <- c( 58.16354, 62.34342, 31.43701, 62.14333,
#                             51.98859, 29.95968, 55.48009, 50.68146, 51.75618, 31.13289,
#                             47.53024, 56.44128)
# height.trait$nonspec<- c( 47.93985, 57.09998, 43.06760, 51.58106,
#                           44.52435, 40.85160, 50.85945, 44.48371, 43.20859, 43.92655, 45.15222, 47.83641)
# # calculate decomposition of total trait variance 
# x1<-trait.flex.anova(~MOWING, specific, nonspec, data=height.trait)
# 
# print(x1)
# par(ask=TRUE)
# plot(x1)
# # calculate decomposition of the factorial model 
# x2<-trait.flex.anova(~MOWING*FERTIL, specific, nonspec, data=height.trait)
# print(x2)
# plot(x2)
# plot(x2, legend.pos='none')
# plot(x2, plot.total=TRUE)
# plot(x2, plot.covar=TRUE)
# # and now the relative values ...
# plot(x1, use.percentage=T)
# plot(x2, use.percentage=T)
# plot(x2, legend.pos=”none”, use.percentage=T)
# plot(x2, plot.total=TRUE, use.percentage=T)
# plot(x2, plot.covar=TRUE, use.percentage=T, main=”Height”) par(ask=FALSE)
# 
# 
# 
# trait.flex.anova(formula, specif.avg, const.avg, ...) 
# 
# 
# x2<-trait.flex.anova(~MOWING*FERTIL, specific, nonspec, data=height.trait)
# x2<-trait.flex.anova(~fish, com_specific_length, overall_length, 
#                      data = count_length_processed2)
# 
# 
# 
# count_length_processed2 %>% 
#   ggplot() +
#   geom_violin(aes(x = fish, y = com_specific_length))
# 
# 
# count_length_processed2 %>% 
#   ggplot() +
#   geom_violin(aes(x = fish, y = overall_length))
# 
# 
# 
# left_join(length_processed, processed_lake_info, by = 'lake') %>% 
#   group_by(fish, taxa) %>% 
#   summarize(sample_size = n()) %>% 
#   pivot_wider(names_from = 'fish', values_from = 'sample_size')
# 
# left_join(length_processed, processed_lake_info, by = 'lake') %>% 
#   group_by(taxa) %>% 
#   summarize(sample_size = n())
# 
# length_mean_info_fish <- left_join(length_processed, processed_lake_info, by = 'lake') %>% 
#   mutate(fish = as.character(fish)) %>% 
#   group_by(taxa) %>% 
#   mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>%
#   ungroup() %>% 
#   group_by(fish, taxa) %>% 
#   summarize(fish_mean_length = mean(length_mm, na.rm = TRUE),
#             overall_mean_length = first(overall_mean_length),
#             .groups = 'drop') %>% 
#   mutate(taxa_fish = paste(taxa, fish, sep = "_"))
# 
# 
# processed_mean_length_fish <- left_join(count_dat_processed1_subset, processed_lake_info, by = 'lake') %>% 
#   mutate(taxa_fish = paste(taxa, fish, sep = "_")) %>%
#   left_join(., length_mean_info_fish, by = 'taxa_fish') %>% 
#   group_by(lake) %>% 
#   summarize(com_specific_length = weighted.mean(fish_mean_length, w = rel_abund),
#             overall_length = weighted.mean(overall_mean_length, w = rel_abund),
#             fish = as.factor(first(fish.x)),
#             .groups = 'drop')
# 
# trait.flex.anova(~fish, com_specific_length, overall_length, 
#                  data = processed_mean_length_fish)
# 
# length_processed %>% 
#   group_by(taxa) %>% 
#   mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   group_by(taxa_lake) %>% 
#   summarize(com_mean_length = mean(length_mm, na.rm = TRUE),
#             overall_mean_length = first(overall_mean_length),
#             .groups = 'drop')
# 
# 
# left_join(count_dat_processed1_subset, 
#           length_processed %>% 
#             group_by(taxa) %>% 
#             mutate(overall_mean_length = mean(length_mm, na.rm = TRUE)) %>% 
#             ungroup() %>% 
#             group_by(taxa_lake) %>% 
#             summarize(com_mean_length = mean(length_mm, na.rm = TRUE),
#                       overall_mean_length = first(overall_mean_length),
#                       .groups = 'drop'), 
#           by = 'taxa_lake') %>%
#   filter(!is.nan(overall_mean_length) | !is.nan(com_mean_length)) %>%
#   group_by(lake) %>% 
#   summarize(com_specific_length = weighted.mean(com_mean_length, w = rel_abund),
#             overall_length = weighted.mean(overall_mean_length, w = rel_abund), 
#             .groups = 'drop')
# 
# subset_length_year
# 
# count_dat_processed1_subset


#ggsave(plot = species_level_shift_plot,
#       filename = here('figures', 'intraspecific_results', 'species_level_shift_plot.pdf'), 
#       width = 10*0.7, height = 11*0.7, device = 'pdf')


### ROUGH EXAMPLE OF HOW TO MODIFY ABOVE PLOT TO ALLOW FOR FACETTING ACROSS MULTIPLE TRAITS
# prac_df <- length_mean_info_fish_init %>%
#   filter(!is.na(length_mm)) %>% 
#   group_by(taxa) %>% 
#   filter(all(c('1', '0') %in% fish)) %>%
#   ungroup() %>% 
#   mutate(fish_info = factor(if_else(fish == '0', 'Fishless', 'Fish'), levels = c('Fishless', 'Fish')),
#          taxa = factor(taxa, levels = taxa_descending_order))
# 
# 
# rbind(prac_df %>% mutate(trait = 'length1'),
#       prac_df %>% mutate(trait = 'length2')) %>% 
#   ggplot(aes(y = taxa, x = length_mm, fill = fish_info, color = fish_info)) + #height = ..count..
#   geom_density_ridges(stat = "density_ridges",
#                       scale = 0.98, rel_min_height = .01, alpha = 0.5,
#                       point_size = 0.5, size = 0.25, #point_shape = "|",
#                       position = position_points_jitter(width = 0.03, height = 0.05, yoffset = -0.03) ) +
#   geom_jitter(data = . %>% 
#                 filter(fish_info == 'Fishless'),
#               #aes(y = as.numeric(taxa) - 0.03),
#               aes(y = as.numeric(taxa) - 0.03),
#               size = 0.3, alpha = 0.7,
#               position = position_jitterdodge(dodge.width = 0,
#                                               jitter.width = 0.07,
#                                               jitter.height = 0)) +
#   geom_jitter(data = . %>% filter(fish_info == 'Fish'),
#               aes(y = as.numeric(taxa) - 0.08),
#               size = 0.3, alpha = 0.7,
#               position = position_jitterdodge(dodge.width = 0,
#                                               jitter.width = 0.07,
#                                               jitter.height = 0)) +
#   geom_point(data = rbind(length_info_processed %>% mutate(trait = 'length1'),
#                           length_info_processed %>% mutate(value = value/5, trait = 'length2')) %>%
#                filter(mean_type == 'Fixed mean'),
#              mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#495057', fill = '#495057') +
#   #geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fishless)'),
#   #           mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#c1121f', fill = '#c1121f') +
#   #geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fish)'),
#   #           mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#662e9b', fill = '#662e9b') +
#   geom_textbox(
#     data = rbind(mannu_results %>% mutate(trait = 'length1'),
#                  mannu_results %>% mutate(trait = 'length2', results_text = paste(results_text, 'green') )),
#     aes(y = taxa, x = Inf, label = results_text),
#     hjust = 1, halign = 0, size = 3,
#     box.colour = NA, fill = NA, # Hide the box fill and outline
#     box.padding = unit(rep(2.75, 4), "pt"), colour = "black",
#     vjust = 1, nudge_y = 0.37, width = NULL
#   )  +
#   #annotate("text", x = 0, y = as.numeric(prac_df$taxa) + 0.5, label = prac_df$test) +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         legend.title = element_blank(),
#         axis.title.y = element_blank()) +
#   xlab('Length (mm)') +
#   facet_wrap(~trait)



# length_mean_info_fish_init %>%
#   filter(!is.na(length_mm)) %>% 
#   group_by(taxa) %>% 
#   filter(all(c('1', '0') %in% fish)) %>% 
#   mutate(fish_info = factor(if_else(fish == '0', 'Fishless', 'Fish'), levels = c('Fishless', 'Fish') )) %>% 
#   ggplot(aes(y = taxa, x = length_mm, fill = fish_info, color = fish_info)) +
#   geom_jitter(size = 0.3, alpha = 0.2, position = position_jitterdodge(jitter.width = 0.7)) +
#   geom_violin(alpha = 0.6, position = position_dodge(width = 0.7), color = 'gray') +
#   theme_bw()

# species_level_shift_plot <- length_mean_info_fish_init %>%
#   filter(!is.na(length_mm)) %>% 
#   group_by(taxa) %>% 
#   filter(all(c('1', '0') %in% fish)) %>%
#   ungroup() %>% 
#   mutate(fish_info = factor(if_else(fish == '0', 'Fishless', 'Fish'), levels = c('Fishless', 'Fish')),
#          taxa = factor(taxa, levels = taxa_descending_order)) %>% 
#   ggplot(aes(y = taxa, x = length_mm, fill = fish_info, color = fish_info)) + #height = ..count..
#   geom_density_ridges(stat = "density_ridges",
#                       scale = 0.98, rel_min_height = .01, alpha = 0.5,
#                       point_size = 0.5, size = 0.25, #point_shape = "|",
#                       position = position_points_jitter(width = 0.03, height = 0.05, yoffset = -0.03) ) +
#   geom_jitter(data = . %>% 
#                 filter(fish_info == 'Fishless'),
#               #aes(y = as.numeric(taxa) - 0.03),
#               aes(y = as.numeric(taxa) - 0.03),
#               size = 0.3, alpha = 0.7,
#               position = position_jitterdodge(dodge.width = 0,
#                                               jitter.width = 0.07,
#                                               jitter.height = 0)) +
#   geom_jitter(data = . %>% filter(fish_info == 'Fish'),
#               aes(y = as.numeric(taxa) - 0.08),
#               size = 0.3, alpha = 0.7,
#               position = position_jitterdodge(dodge.width = 0,
#                                               jitter.width = 0.07,
#                                               jitter.height = 0)) +
#   geom_point(data = length_info_processed %>% filter(mean_type == 'Fixed mean'),
#              mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#495057', fill = '#495057') +
#   geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fishless)'),
#              mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#c1121f', fill = '#c1121f') +
#   geom_point(data = length_info_processed %>% filter(mean_type == 'Mean (fish)'),
#              mapping = aes(y = as.numeric(taxa) + 0.02, x = value), size = 2.25, shape = 1, stroke = 1.4, color = '#662e9b', fill = '#662e9b') +
#   geom_textbox(
#     data = mannu_results,
#     aes(y = taxa, x = Inf, label = results_text),
#     hjust = 1, halign = 0, size = 3,
#     box.colour = NA, fill = NA, # Hide the box fill and outline
#     box.padding = unit(rep(2.75, 4), "pt"), colour = "black",
#     vjust = 1, nudge_y = 0.37, width = NULL
#   )  +
#   #annotate("text", x = 0, y = as.numeric(prac_df$taxa) + 0.5, label = prac_df$test) +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         legend.title = element_blank(),
#         axis.title.y = element_blank()) +
#   xlab('Length (mm)')

#ggsave(plot = vardecomp_fish_vs_fishless_plot,
#       filename = here('figures', 'intraspecific_results', 'vardecomp_fish_vs_fishless_plot.pdf'), 
#       width = 15*0.7, height = 8*0.7, device = 'pdf')


# vardecomp_fish_vs_fishless %>% 
#   filter(quantity %in% c('between_var_prop', 'within_var_prop')) %>% 
#   ggplot() +
#   geom_bar(aes(x = fish_info, y = observed, fill = quantity), stat = 'identity', position = position_dodge()) +
#   geom_errorbar(aes(x = fish_info, ymin = lower_ci, ymax = upper_ci, group = quantity),
#                 width = 0.2, position = position_dodge()) +
#   ylim(0, 1) +
#   theme_bw()

