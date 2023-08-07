####################################
####################################
### ANLAYSIS OF ZOOPLANKTON DATA ###
####################################
####################################

#https://github.com/ShanKothari/DecomposingFD
#https://www.biorxiv.org/content/10.1101/530782v1

#COUNT DATA QUESTIONS
# #keratella, keratella 2 spine, keratella 1 spine in 2018
# #should d. middenorfiana/pulex and d. middenorfiana be treated separately or be lumped together?
# #difference between immature and juveline?
# 
# #Cyclapoid vs Cyclopoid unident. --> from the same lake


##########################
### SCRIPT PREPARATION ###
##########################

### load libraries ###
library(here)
library(tidyverse)
library(gawdis)
library(readxl)
library(picante) #matrix randomization

#load custom functions
source(here("scripts", "custom_functions_zooplankton_project.R"))


### LOADING DATA ###
#(processed) count and length data
count_dat_processed_list <- readRDS(file = here('data', 
                                                'processed_data', 
                                                'count_dat_processed_list.rds'))

zoop_length_total_list0 <- readRDS(file = here('data', 
                                               'processed_data', 
                                               'zoop_length_total_list0.rds'))

#fish presence/absence info
lindsey_master_data <- as.data.frame(read_excel(here('data', 
                                                     'WR_data_ALL_Boyle.xlsx')))

#trait data (not including length data)
trait_info_init <- as.data.frame(read_excel(here('data', 
                                                 'Zooplankton_eggs_and_length_5_9_2022.xlsx'), 
                                            sheet = 'Inter Traits'))

all_taxa_traits <- readRDS(here('data', 
                                'processed_data', 
                                'all_taxa_trait_info_processed.rds'))



##########################
### PROCESSING OF DATA ###
##########################
#*** most processing done in zooplankton_data_processing_script.R

###  FISH PRESENCE/ABSENCE INFO FOR EACH LAKE ###
fish_info <- lindsey_master_data %>% 
  select(lake, fish)

processed_lake_info <- fish_info[!duplicated(fish_info$lake),] %>%
  mutate(lake = recode(lake,
                       'lake37' = "lake_37",
                       'blackjoe' = "black_joe",
                       'upper_blackjoe' = "upper_black_joe",
                       "below_tb" = "lake_below_tb")) #%>% 
  #select(!lake)


#remove euchlanis and filinia and move secies to a row name
all_taxa_traits_processed_init <- all_taxa_traits %>% 
  filter(!(species %in% c('euchlanis', 'filinia')) ) %>%
  column_to_rownames(var = 'species')

count_dat_processed_init <- count_dat_processed_list %>%
  bind_rows() %>% 
  filter(Lake != "Lost Lake") %>% #remove lost lake
  mutate(Lake = dplyr::recode(Lake, "N of Blue" = "North Blue")) %>%
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



#####################################################################
### COMMUNITY-LEVEL MULTIDIM FD ANALYSIS (FRic, FDis, FDiv, FEve) ###
#####################################################################

### STEP 1: CALCULATE DISSIMILARITY MATRIX ###
#de Bello et al.'s (2021) extension of Gower distance used for dissim. mat. construction
zoo_dist <- dist_wrapper(df = all_taxa_traits_processed,
                         method = c('dist', 'gowdis', 'gawdis')[3], 
                         arg_list = list(w.type = 'optimized', opti.maxiter = 500)) 

### STEP 2: DETERMINE OPTIMAL FUNCTIONAL SPACE ###
zoo_msd_calc <- eval_func_space(trait_data = zoo_dist,
                                   dim = 10,
                                   eval_option = c("corank", "mad"),
                                   show_messages = TRUE,
                                   rand_rep = 99)

#information about the optimal functional space
### *** WITH THE FINAL DATASET, THIS SHOULD NOT BE TRUNCATED AT 2 DIMENSIONS
optimal_func_space_info <- zoo_msd_calc$func_eval$corank[which(zoo_msd_calc$func_eval$corank$dim <= 4),] %>% 
  filter(auc == max(auc))

#process community data based on optimal functional space
count_dat_processed <- count_dat_processed_init %>%
  #filter(taxa %in% focal_taxa) %>%
  group_by(lake) %>% 
  filter(n() >= (optimal_func_space_info$dim + 1) ) %>% 
  ungroup() %>%
  pivot_wider(names_from = lake, values_from = total_abund, values_fill = 0) %>% 
  as.data.frame()


#check to make sure that the processed community dataset still has the identical species as the trait dataset
if (identical(sort(count_dat_processed$taxa), sort(rownames(all_taxa_traits_processed)))) {
  message('The processed community dataset (removing communities with too few species) has the same species as the trait dataset. Analysis can proceed!')
} else {
  stop('The processed community dataset (removing communities with too few species) does NOT have the same species as the trait dataset. Please address before proceeding')
}


### STEP 3: CALCULATE OBSERVED FUNCTIONAL DIVERSITY ###
fdcontr_calc <- calc_fd_contr_inputdist(data = count_dat_processed,
                                        trait_dist = zoo_dist,
                                        species = 'all',
                                        metrics = c("FRic", "FDis", "FDiv", "FEve"),
                                        std_traits = TRUE,
                                        pcoa_dims = optimal_func_space_info$dim,
                                        correction = optimal_func_space_info$correction,
                                        feve_calc_type = c("reduced_PCoA", "FD")[1],
                                        mst_package = c('vegan', 'ape')[1],
                                        w_abund = TRUE,
                                        calc_contr = FALSE,
                                        show_messages = TRUE,
                                        include_absent_species = FALSE,
                                        unique_id = TRUE,
                                        output_PCoA = TRUE)


### STEP 4: CALCULATE STANDARDIZED EFFECT SIZE OF THE METRICS ###

#create 1000 randomized communities
set.seed(482938235)
randomized_com_list <- lapply(seq(1000), function(x) {
  randomizeMatrix(t(count_dat_processed %>% column_to_rownames('taxa')), 
                  null.model = c("frequency", "richness", "independentswap", "trialswap")[3], 
                  iterations = 1000) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'taxa')
})

#calculate the FD metrics on each randomized dataset
randomized_fd_calc_progress <- txtProgressBar(min = 0, max = length(randomized_com_list), 
                                              char = "*", style = 3)

rand_calcs_list <- list()
for (i in seq_len(length(randomized_com_list))) {
  rand_calcs_list[[i]] <- calc_fd_contr_inputdist(data = randomized_com_list[[i]],
                                                  trait_dist = zoo_dist,
                                                  species = 'all',
                                                  metrics = c("FRic", "FDis", "FDiv", "FEve"),
                                                  std_traits = TRUE,
                                                  pcoa_dims = optimal_func_space_info$dim,
                                                  correction = optimal_func_space_info$correction,
                                                  feve_calc_type = c("reduced_PCoA", "FD")[1],
                                                  mst_package = c('vegan', 'ape')[1],
                                                  w_abund = TRUE,
                                                  calc_contr = FALSE,
                                                  show_messages = FALSE,
                                                  include_absent_species = FALSE,
                                                  unique_id = TRUE,
                                                  output_PCoA = TRUE)$FD_dataframe
  
  setTxtProgressBar(randomized_fd_calc_progress, i)
}

## also easy to do with lapply but not as easy to monitor progress w/ a progress bar
# rand_calcs_list <- lapply(randomized_com_list, function(RANDOM_COM) {
#   
#   calc_fd_contr_inputdist(data = RANDOM_COM,
#                           trait_dist = zoo_dist,
#                           species = 'all',
#                           metrics = c("FRic", "FDis", "FDiv", "FEve"),
#                           std_traits = TRUE,
#                           pcoa_dims = optimal_func_space_info$dim,
#                           correction = optimal_func_space_info$correction,
#                           feve_calc_type = c("reduced_PCoA", "FD")[1],
#                           mst_package = c('vegan', 'ape')[1],
#                           w_abund = TRUE,
#                           calc_contr = FALSE,
#                           show_messages = FALSE,
#                           include_absent_species = FALSE,
#                           unique_id = TRUE,
#                           output_PCoA = TRUE)$FD_dataframe
#   }
# )

#calculate SES for each metric
fd_ses_df <- rand_calcs_list %>% 
  bind_rows() %>% 
  group_by(community) %>% 
  summarize(fric_mean = mean(full_FRic),
         fric_sd = sd(full_FRic),
         fdis_mean = mean(full_FDis),
         fdis_sd = sd(full_FDis),
         fdiv_mean = mean(full_FDiv),
         fdiv_sd = sd(full_FDiv),
         feve_mean = mean(full_FEve),
         feve_sd = sd(full_FEve)) %>% 
  ungroup() %>% 
  select(community, fric_mean, fric_sd, fdis_mean , fdis_sd, fdiv_mean, fdiv_sd, feve_mean, feve_sd) %>% 
  left_join(., fdcontr_calc$FD_dataframe, by = 'community') %>% 
  mutate(fric_ses = (full_FRic - fric_mean)/fric_sd,
         fdis_ses = (full_FDis - fdis_mean)/fdis_sd,
         fdiv_ses = (full_FDiv - fdiv_mean)/fdiv_sd,
         feve_ses = (full_FEve - feve_mean)/feve_sd) %>% 
  select(community, fric_ses, fdis_ses, fdiv_ses, feve_ses) %>% 
  rename(lake = community) #%>%
  # mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
  # mutate(lake = recode(lake,
  # 'annika_lake' = 'annika',
  # "bear_west" = "west_bear",
  #  "east_no_name" = "east_noname",
  #  "frozen_lake_1" = "frozen1",
  #  "frozen_lake_2" = "frozen2",
  #  "island_lake" = 'island',
  #  "lindsey_lake" = "lindsey",
  #  "little_mountain_sheep" = "little_mtnsheep",
  #  "low_deep_creek" = "lower_deepcreek",
  #  "lower_indian_basin" = "lower_indian",
  #  "mid_deep_creek" = "mid_deepcreek",
  #  "mid_indian_basin" = "mid_indian",
  #  "n_of_blue" = "north_blue",
  #  "no_name_west" = "west_noname",
  #  "upper_indian_basin" = "upper_indian",
  #  "west_of_mt_lester" = "west_mtlester"))

#n of blue switched to north blue but north blue was already in the dataset


write.csv(fd_ses_df, file = here('results', 'alpha_ses_fd.csv'), row.names = FALSE)




fd_ses_df_lake <- left_join(fd_ses_df, processed_lake_info, by = 'lake')  %>%
  #filter(lake != 'lost_lake') %>% 
  #select(!community) %>% 
  pivot_longer(cols = ends_with('ses'),
               names_to = "fd_metric", 
               values_to = "value") %>% 
  mutate(fish = as.character(fish)) %>% 
  



#######################################
### VISUALIZING MULTIDIM FD METRICS ###
#######################################

#library(ggsignif)

fd_ses_df_lake_plot <- fd_ses_df_lake %>%
  mutate(fish = if_else(fish == '1', 'fish', 'fishless'),
         metric = factor(recode(fd_metric,
                                'fric_ses' = 'Richness',
                                'fdis_ses' = 'Dispersion',
                                'fdiv_ses' = 'Divergence',
                                'feve_ses' = 'Evenness'),
                         levels = c('Richness', 'Dispersion', 'Divergence', 'Evenness')),
         fish = factor(fish, levels = c('fishless', 'fish'))) %>%
  ggplot() +
  #geom_boxplot(aes(x = metric, y = value, fill = fish),
  #             alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
  geom_violin(aes(x = metric, y = value, fill = fish),
              alpha = 0.3, color = '#a6a6a6', lwd = 0.8) +
  geom_point(aes(x = metric, y = value, fill = fish, group = fish, color = fish),
             size = 2,
             position = position_jitterdodge(dodge.width = 0.875, jitter.width = 0.2)) +
  xlab('Functional diversity') +
  ylab('Standardized effect size') +
  scale_fill_manual(values = c('#DC3220', '#005AB5')) +
  scale_color_manual(values = c('#DC3220', '#005AB5')) +
  theme_bw() +
  theme(axis.title = element_text(size = 19),
        axis.text.x = element_text(size = 15),
        legend.title=element_blank(),
        legend.text = element_text(size = 15)) 

ggsave(filename = here('figures', 'updated_multidim_FD_violin.png'),
       plot = fd_ses_df_lake_plot,
       width = 13*0.8, height = 7*0.8)



###############################################
### SCHEINER ET AL. (2016) ALPHA FD METRICS ###
###############################################

alpha_ftd <- FTD.comm(tdmat = zoo_dist, 
                      spmat = count_dat_processed %>% 
                        column_to_rownames(var = 'taxa') %>% 
                        t(), 
                     q = 1, abund = TRUE, 
                     match.names = TRUE)$com.FTD %>% 
  rename(lake = community) #%>%
  # mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
  # mutate(lake = recode(lake, 
  #                      'annika_lake' = 'annika',
  #                      "bear_west" = "west_bear",
  #                      "east_no_name" = "east_noname",
  #                      "frozen_lake_1" = "frozen1",
  #                      "frozen_lake_2" = "frozen2",
  #                      "island_lake" = 'island',
  #                      "lindsey_lake" = "lindsey",
  #                      "little_mountain_sheep" = "little_mtnsheep",
  #                      "low_deep_creek" = "lower_deepcreek",
  #                      "lower_indian_basin" = "lower_indian",
  #                      "mid_deep_creek" = "mid_deepcreek",
  #                      "mid_indian_basin" = "mid_indian",
  #                      "n_of_blue" = "north_blue",
  #                      "no_name_west" = "west_noname",
  #                      "upper_indian_basin" = "upper_indian",
  #                      "west_of_mt_lester" = "west_mtlester"))

write.csv(alpha_ftd, file = here('results', 'alpha_scheiner_metrics.csv'), row.names = FALSE)


processed_scheiner_data_vis <- left_join(alpha_ftd, 
                                         processed_lake_info, 
                                         by = 'lake')  %>%
  #select(!ccommunity) %>% 
  filter(lake != 'lost_lake') %>% 
  #select(!c(community, q) ) %>% 
  pivot_longer(cols = !c(lake, fish),
               names_to = "fd_metric", 
               values_to = "value") %>% 
  mutate(fish_info = factor(if_else(fish == 1, 'fish', 'fishless'), 
                            levels = c('fishless', 'fish')) )  %>% 
  #mutate(fish = as.character(fish)) %>% 
  filter((fd_metric %in% c('nsp', 'M.prime', 'qEt', 'qDTM') ) ) %>%
  mutate(fd_metric = factor(fd_metric, levels = c('nsp', 'M.prime', 'qEt', 'qDTM')),
         fd_metric = recode(fd_metric, 'nsp' = 'Species richness') )

ggplot(data = processed_scheiner_data_vis) +
  geom_boxplot(aes(x = fd_metric, y = value, fill = fish_info),
               alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
  geom_point(aes(x = fd_metric, y = value, fill = fish_info, group = fish_info, color = fish_info),
             size = 2,
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
  #geom_histogram(aes(x = value, fill = fish)) +
  facet_wrap(vars(fd_metric), scales = "free") +
  theme_bw() +
  scale_fill_manual(values = c('#DC3220', '#005AB5')) +
  scale_color_manual(values = c('#DC3220', '#005AB5')) +
  theme(legend.title = element_blank()) +
  xlab('Metric') + ylab('Value')

ggsave(filename = here('figures', 'scheiner_fd_metrics_boxplot.png'), 
       width = 13*0.8, height = 9*0.8)



##################################
### SUMMARY TABLE OF FD VALUES ###
##################################

fd_summary_table <- left_join(alpha_ftd, 
          processed_lake_info, 
          by = 'lake')  %>%
  left_join(., 
            fd_ses_df, 
            by = 'lake') %>%
  left_join(., 
            fdcontr_calc$FD_dataframe %>% 
              rename(lake = community) %>% 
              select(!species_richness), 
            by = 'lake') %>%
  rename_with(~ gsub("_", " ", .x)) %>% 
  group_by(fish) %>% 
  select(-c(lake, q)) %>% 
  summarise(across(everything(), 
                   list(min = min, 
                        max = max,
                        mean = mean,
                        sd = sd)
  )
  ) %>% 
  pivot_longer(cols = !fish, 
               names_to = c('metric', 'summary'),
               names_sep = "_",
               values_to = "value") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(fish = if_else(fish == 0, 'no', 'yes')) %>% 
  rename_with(~ gsub(" ", "_", .x))

write.csv(fd_summary_table, 
          file = here('results', 'summary_info', 'fd_summary_table.csv'), 
          row.names = FALSE)



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# prac_ftd_processed <- prac_ftd %>%
#   mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
#   mutate(lake = recode(lake, 
#                        'annika_lake' = 'annika',
#                        "bear_west" = "west_bear",
#                        "east_no_name" = "east_noname",
#                        "frozen_lake_1" = "frozen1",
#                        "frozen_lake_2" = "frozen2",
#                        "island_lake" = 'island',
#                        "lindsey_lake" = "lindsey",
#                        "little_mountain_sheep" = "little_mtnsheep",
#                        "low_deep_creek" = "lower_deepcreek",
#                        "lower_indian_basin" = "lower_indian",
#                        "mid_deep_creek" = "mid_deepcreek",
#                        "mid_indian_basin" = "mid_indian",
#                        "n_of_blue" = "north_blue",
#                        "no_name_west" = "west_noname",
#                        "upper_indian_basin" = "upper_indian",
#                        "west_of_mt_lester" = "west_mtlester"))

# left_join(alpha_ftd %>% rename(lake = community), processed_lake_info, by = 'lake')  %>%
#   #select(!ccommunity) %>% 
#   filter(lake != 'lost_lake') %>% 
#   #select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish)) %>% 
#   #filter(fd_metric %in% c('qDTM', 'qEt')) %>%
#   filter(fd_metric %in% c('M.prime')) %>%
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
#   theme_bw()
# 
# 
# left_join(alpha_ftd %>% rename(lake = community), processed_lake_info, by = 'lake') %>%
#   #select(!ccommunity) %>% 
#   filter(lake != 'lost_lake') %>% 
#   #select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish)) %>% 
#   #filter(fd_metric %in% c('qDTM', 'qEt')) %>%
#   filter(fd_metric %in% c('M.prime')) %>%
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
#   theme_bw()
# 
# 
# left_join(alpha_ftd %>% rename(lake = community), processed_lake_info, by = 'lake')  %>%
#   #select(!ccommunity) %>% 
#   filter(lake != 'lost_lake') %>% 
#   #select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish)) %>% 
#   #filter(fd_metric %in% c('qDTM', 'qEt')) %>%
#   #filter(fd_metric %in% c('M.prime')) %>%
#   filter(!(fd_metric %in% c('M', 'nsp', 'q') ) ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
#   theme_bw()


# left_join(alpha_ftd %>% rename(lake = community), processed_lake_info, by = 'lake')  %>%
#   #select(!ccommunity) %>% 
#   filter(lake != 'lost_lake') %>% 
#   #select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish)) %>% 
#   #filter(fd_metric %in% c('qDTM', 'qEt')) %>%
#   #filter(fd_metric %in% c('M.prime')) %>%
#   filter(!(fd_metric %in% c('M', 'nsp', 'q') ) ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
#   facet_wrap(vars(fd_metric)) +
#   theme_bw()


# ggsave(filename = here('figures', 'updated_multidim_FD_violin.png'), 
#        width = 13*0.8, height = 7*0.8)



# prac_ftd <- FTD.comm(tdmat = zoo_dist, 
#          spmat = count_dat_processed %>% 
#            column_to_rownames(var = 'taxa') %>% 
#            t(), 
#          q=1, abund = TRUE, 
#          match.names = TRUE)$com.FTD %>%
#   mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
#   mutate(lake = recode(lake, 
#                        'annika_lake' = 'annika',
#                        "bear_west" = "west_bear",
#                        "east_no_name" = "east_noname",
#                        "frozen_lake_1" = "frozen1",
#                        "frozen_lake_2" = "frozen2",
#                        "island_lake" = 'island',
#                        "lindsey_lake" = "lindsey",
#                        "little_mountain_sheep" = "little_mtnsheep",
#                        "low_deep_creek" = "lower_deepcreek",
#                        "lower_indian_basin" = "lower_indian",
#                        "mid_deep_creek" = "mid_deepcreek",
#                        "mid_indian_basin" = "mid_indian",
#                        "n_of_blue" = "north_blue",
#                        "no_name_west" = "west_noname",
#                        "upper_indian_basin" = "upper_indian",
#                        "west_of_mt_lester" = "west_mtlester"))
# 
# prac_ftd_processed <- prac_ftd %>%
#   mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
#   mutate(lake = recode(lake, 
#                        'annika_lake' = 'annika',
#                        "bear_west" = "west_bear",
#                        "east_no_name" = "east_noname",
#                        "frozen_lake_1" = "frozen1",
#                        "frozen_lake_2" = "frozen2",
#                        "island_lake" = 'island',
#                        "lindsey_lake" = "lindsey",
#                        "little_mountain_sheep" = "little_mtnsheep",
#                        "low_deep_creek" = "lower_deepcreek",
#                        "lower_indian_basin" = "lower_indian",
#                        "mid_deep_creek" = "mid_deepcreek",
#                        "mid_indian_basin" = "mid_indian",
#                        "n_of_blue" = "north_blue",
#                        "no_name_west" = "west_noname",
#                        "upper_indian_basin" = "upper_indian",
#                        "west_of_mt_lester" = "west_mtlester"))
# 
# left_join(prac_ftd_processed, processed_lake_info, by = 'lake')  %>%
#   #select(!ccommunity) %>% 
#   filter(lake != 'lost_lake') %>% 
#   select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish)) %>% 
#   #filter(fd_metric %in% c('qDTM', 'qEt')) %>%
#   filter(fd_metric %in% c('M.prime')) %>%
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0)) +
#   theme_bw()
# 
# 
# prac_ftd_processed %>% 
#   filter(lake != 'lost_lake') %>% 
#   left_join(., processed_lake_info, by = 'lake') %>% 
#   mutate(fish = if_else(fish == '1', 'fish', 'fishless')  ) %>% 
#   select(qEt, qDTM, fish) %>% 
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0))
# 
#   


#M --> a measure of the magnitude of dispersion
#      it has a range of [0,(S-1)/S] and measures standardized distance.

#qHT --> measure of variability among pairwise distances with the Hill function

#qDT --> indicates the effective number of equally distant species and has a range of [1,S]
#        qDT can be decomposed into components of richness and evenness:

#qDTM --> functional-trait species diversity andmean dispersion into an integrated metric, functional-trait dispersion
#         the metric combines species richness, functional-trait evennessand mean dispersion





# prac_ftd_metdat <- left_join(prac_ftd, fish_info[!duplicated(fish_info$lake),], by = 'lake')  %>%
#   #select(!ccommunity) %>% 
#   select(!c(community, q) ) %>% 
#   pivot_longer(cols = !c(lake, fish),
#                names_to = "fd_metric", 
#                values_to = "value") %>% 
#   mutate(fish = as.character(fish))
# 
# prac_ftd_metdat 
# 
# 
# 
# prac_ftd_metdat %>% 
#   mutate(fish = if_else(fish == '1', 'fish', 'fishless')  ) %>% 
#   #filter(fd_metric %in% c('M.prime', 'qDTM', 'qEt')) %>% 
#   filter(fd_metric %in% 'qEt') %>% 
#   ggplot() +
#   geom_boxplot(aes(x = fd_metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = fd_metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0))
# 
# 
# fd_ses_df_lake %>% 
#   mutate(fish = if_else(fish == '1', 'fish', 'fishless'),
#          metric = factor(recode(fd_metric,
#                          'fric_ses' = 'Richness',
#                          'fdis_ses' = 'Dispersion',
#                          'fdiv_ses' = 'Divergence',
#                          'feve_ses' = 'Evenness'),
#                          levels = c('Richness', 'Dispersion', 'Divergence', 'Evenness')),
#          fish = factor(fish, levels = c('fishless', 'fish'))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = metric, y = value, fill = fish),
#                alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2))



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# #THIS CAUSES MY R SESSION TO CRASH
# library(betapart)
# 
# count_dat_processed_beta_1 <- count_dat_processed %>%
#   column_to_rownames(var = 'taxa') %>%
#   t() %>%
#   as.data.frame()
# 
# count_dat_processed_beta <- count_dat_processed_beta_1 %>% 
#   filter(!(rownames(count_dat_processed_beta_1) %in% c('East Bear', 
#                                                        'Big Stough', 
#                                                        'Black Joe', 
#                                                        'Frozen Lake 1',
#                                                        'Bluebell',
#                                                        'Boulder',
#                                                        'Casey')) )
# 
# #count_dat_processed_beta_1[rownames(count_dat_processed_beta_1) == 'Little Island',]
# 
# count_dat_processed_beta[count_dat_processed_beta != 0] <- 1
# count_dat_processed_beta_final <- count_dat_processed_beta[,match(rownames(zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4]), colnames(count_dat_processed_beta))]
# 
# pairwise_lake_combos <- as.data.frame(t(combn(rownames(count_dat_processed_beta_final), 2)))
# colnames(pairwise_lake_combos) <- c('lake1', 'lake2')
# 
# pairwise_lake_combos_step2 <- pairwise_lake_combos %>% 
#   rowwise() %>% 
#   mutate(lake_combo_id = paste(sort(c(lake1, lake2)), collapse = "_"))
# 
# lake_combo_remove_list <- list(c('East Bear', 'Rapid'),
#                                c('47Aaron', 'Elephant Head'),
#                                c('Bear West', 'Boulder'),
#                                c("Big Stough", 'Bluebell'),
#                                c("Big Stough", 'East Bear'),
#                                c("Big Stough", 'Frozen Lake 1'),
#                                c("Big Stough", 'Island Lake'),
#                                c('Black Joe', 'Elephant Head'),
#                                c('Bluebell', 'Casey'),
#                                c('Bluebell', 'Btwn Wrbnt Bldr'),
#                                c('Boulder', 'Btwn Wrbnt Bldr'),
#                                c('Boulder', 'Casey'),
#                                c('Boot', 'Windy'),
#                                c('Boulder', 'Crater'),
#                                c('Casey', 'Little Island'),
#                                c('Casey', 'Island Lake'),
#                                c('Cliff', 'Elk'),
#                                c('Cliff', 'Dome'))
# 
# 
# remove_combo_vec <- sapply(lake_combo_remove_list, function(x) paste0(sort(c(x[1], x[2])), collapse = "_")  )
# 
# pairwise_lake_combos_step3 <- pairwise_lake_combos_step2 %>% 
#   filter(!(lake_combo_id %in% remove_combo_vec) )
# 
# 
# beta_list <- list()
# for (i in seq_len(nrow(pairwise_lake_combos_step3))) {
#   lake1 <- pairwise_lake_combos_step3[i,'lake1']
#   lake2 <- pairwise_lake_combos_step3[i,'lake2']
#   message('index: ', i, '; lake1: ', lake1, '; lake2: ', lake2)
#   
#   count_dat_processed_beta_step1 <- count_dat_processed %>%
#     column_to_rownames(var = 'taxa') %>%
#     t() %>%
#     as.data.frame()
#   
#   count_dat_processed_beta_subset <- count_dat_processed_beta_step1 %>% 
#     filter(row.names(count_dat_processed_beta_step1) %in% c(lake1, lake2))
#     
#   count_dat_processed_beta_subset[count_dat_processed_beta_subset != 0] <- 1
#   count_dat_processed_beta_step2 <- count_dat_processed_beta_subset[,match(rownames(zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4]), colnames(count_dat_processed_beta_subset))]
#   
#   beta_list[[paste(lake1, lake2, sep = '_')]] <- functional.beta.pair(x = count_dat_processed_beta_step2,
#                        traits = zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4],
#                        index.family = 'sorensen')
#   
# }

# count_dat_processed_beta_step1 <- count_dat_processed %>%
#   column_to_rownames(var = 'taxa') %>%
#   t() %>%
#   as.data.frame()
# 
# count_dat_processed_beta_subset <- count_dat_processed_beta_step1 %>% 
#   filter(row.names(count_dat_processed_beta_step1) %in% c("Big Stough", 'Frozen Lake 1'))
# 
# count_dat_processed_beta_subset[count_dat_processed_beta_subset != 0] <- 1
# count_dat_processed_beta_step2 <- count_dat_processed_beta_subset[,match(rownames(zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4]), colnames(count_dat_processed_beta_subset))]
# 
# functional.beta.pair(x = count_dat_processed_beta_step2,
#                                                                     traits = zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4],
#                                                                     index.family = "sorensen")
# 
# 
# 
# functional.beta.pair(x = count_dat_processed_beta_final,
#                      traits = zoo_msd_calc$supplementary_info$PCoA_list$quasieuclid[,1:4],
#                      index.family = 'sorensen')




###############################
### JASM POSTER MAIN RESULT ###
###############################

# library(ggsignif)
# 
# fd_ses_df_lake %>%
#   mutate(fish = if_else(fish == '1', 'fish', 'fishless'),
#          metric = factor(recode(fd_metric,
#                          'fric_ses' = 'Richness',
#                          'fdis_ses' = 'Dispersion',
#                          'fdiv_ses' = 'Divergence',
#                          'feve_ses' = 'Evenness'),
#                          levels = c('Richness', 'Dispersion', 'Divergence', 'Evenness')),
#          fish = factor(fish, levels = c('fishless', 'fish'))) %>%
#   ggplot() +
#   #geom_boxplot(aes(x = metric, y = value, fill = fish),
#   #             alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
#   geom_violin(aes(x = metric, y = value, fill = fish),
#                alpha = 0.3, color = '#a6a6a6', lwd = 0.8) +
#   geom_point(aes(x = metric, y = value, fill = fish, group = fish, color = fish),
#              size = 2,
#              position = position_jitterdodge(dodge.width = 0.875, jitter.width = 0.2)) +
#   xlab('Functional diversity') +
#   ylab('Standardized effect size') +
#   scale_fill_manual(values = c('#DC3220', '#005AB5')) +
#   scale_color_manual(values = c('#DC3220', '#005AB5')) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 19),
#         axis.text.x = element_text(size = 15),
#         legend.title=element_blank(),
#         legend.text = element_text(size = 15)) 
# 
# ggsave(filename = here('figures', 'updated_multidim_FD_violin.png'), 
#        width = 13*0.8, height = 7*0.8)
#  


# prac_ftd <- FTD.comm(tdmat = zoo_dist, 
#                      spmat = count_dat_processed %>% 
#                        column_to_rownames(var = 'taxa') %>% 
#                        t(), 
#                      q = 1, abund = FALSE, 
#                      match.names = TRUE)$com.FTD %>%
#   mutate(lake = gsub(" ", "_", str_trim(tolower(community)))) %>%
#   mutate(lake = recode(lake, 
#                        'annika_lake' = 'annika',
#                        "bear_west" = "west_bear",
#                        "east_no_name" = "east_noname",
#                        "frozen_lake_1" = "frozen1",
#                        "frozen_lake_2" = "frozen2",
#                        "island_lake" = 'island',
#                        "lindsey_lake" = "lindsey",
#                        "little_mountain_sheep" = "little_mtnsheep",
#                        "low_deep_creek" = "lower_deepcreek",
#                        "lower_indian_basin" = "lower_indian",
#                        "mid_deep_creek" = "mid_deepcreek",
#                        "mid_indian_basin" = "mid_indian",
#                        "n_of_blue" = "north_blue",
#                        "no_name_west" = "west_noname",
#                        "upper_indian_basin" = "upper_indian",
#                        "west_of_mt_lester" = "west_mtlester"))

#+
  #geom_signif(stat="identity",
  #           data=data.frame(x=c(0.875, 1.875, 2.875, 3.875),
  #                           xend=c(1.125, 2.125, 3.125, 4.125),
  #                           y=c(2.2, 2.2, 2.2, 2.2),
  #                           annotation = c(" ", " ", " ", " ")),
  #           aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  # geom_text(data = data.frame(x=c(1, 2, 3, 4),
  #                             y=c(2.25, 2.32, 2.41, 2.41),
  #                             annotation = c("***", "NS", ".", ".")),
  #           aes(x = x, y = y, label = annotation), size = c(4, 5, 10, 10) )
  #ylim(-2.5, 2.5)
# 
# ggsave(filename = here('JASM_2022_poster', 'FD_poster_boxplot.png'), 
#        width = 13*0.8, height = 7*0.8)
# 
# summary(aov(value ~ fish, data = fd_ses_df_lake %>% filter(fd_metric == 'fric_ses')))
# summary(aov(value ~ fish, data = fd_ses_df_lake %>% filter(fd_metric == 'fdis_ses')))
# summary(aov(value ~ fish, data = fd_ses_df_lake %>% filter(fd_metric == 'fdiv_ses')))
# summary(aov(value ~ fish, data = fd_ses_df_lake %>% filter(fd_metric == 'feve_ses')))
# 
# 
# summary(lm(value ~ fish, data = fd_ses_df_lake %>% filter(fd_metric == 'fdis_ses')))
# 
# 
# fd_ses_df_lake[!duplicated(fd_ses_df_lake$lake),] %>% 
#   mutate(fish = as.character(fish)) %>% 
#   group_by(fish) %>% 
#   summarize(count = n())
# 
# summary(lm(value ~ fish, data = fd_ses_df_lake[fd_ses_df_lake$fd_metric == 'feve_ses',]))
#   
# table(fd_ses_df_lake$lake, as.character(fd_ses_df_lake$fish))
# 
# 
# sort(unique(fd_ses_df$lake)[!(unique(fd_ses_df$lake) %in% fish_info$lake)])
# sort(unique(fish_info$lake))
# 
# 
# #Hesperodiaptomus shoshone
# 
# egg_len_dat_processed %>% 
#   filter(species == 'H_shoshone') %>% 
#   ggplot() + 
#   geom_histogram(data = egg_len_dat_processed %>% 
#                    filter(species == 'H_shoshone') %>% filter(contains_fish == 'yes'), 
#                  aes(size_mm, fill = contains_fish), 
#                  bins = 10, alpha = 0.4) + 
#   geom_histogram(data = egg_len_dat_processed %>% 
#                    filter(species == 'H_shoshone') %>% filter(contains_fish == 'no'), 
#                  aes(size_mm, fill = contains_fish, alpha = contains_fish), 
#                  bins = 10, alpha = 0.4) +
#   scale_fill_manual(name = NULL, 
#                     values = c('#DC3220', '#005AB5'),
#                     guide = guide_legend(override.aes = list(alpha = 0.4) ) ) +
#   #guides(fill = guide_legend(override.aes = list(alpha = 0.4))) +
#   theme_bw() +
#   theme(legend.position = 'none',
#         axis.title = element_text(size = 18)) +
#   xlab('Size at reproduction (mm)') +
#   ylab('Count')
# 
# ggsave(filename = here('JASM_2022_poster', 'h_shoshone_length_at_repro.png'), 
#        width = 12*0.6, height = 11*0.6)




##############################
### VARIANCE DECOMPOSITION ###
##############################

#[1] "Upper Lost"            "Elephant Head"        
#[3] "Lindsey Lake"          "West of Mt Lester"    
#[5] "Upper Indian Basin"    "Cliff"                
#[7] "Annika Lake"           "Arrowhead"  

# prac_var_decomp <- lapply(setNames(nm = unique(zoop_length_total_list0$length_dat_2018$Lake)), function(LAKE, length_info, abund_info) {
#   
#   length_dat_example <- length_info %>% 
#     filter(Lake == LAKE) %>% 
#     select(taxa, length_mm) %>% 
#     rename(species = taxa,
#            trait = length_mm)
#   
#   total_count <- abund_info %>% 
#     filter(Lake == LAKE) %>% 
#     filter(taxa %in% unique(length_dat_example$species)) %>% 
#     pull(total) %>% 
#     sum()
#   
#   count_dat_example <- abund_info %>% 
#     filter(Lake == LAKE) %>% 
#     filter(taxa %in% unique(length_dat_example$species)) %>% 
#     mutate(rel_abund = total/total_count) %>% 
#     select(taxa, rel_abund) %>% 
#     rename(species = taxa)
#   
#   var_info_df <- var_decomp(abundance_df = count_dat_example, trait_df = length_dat_example)
#   
#   # return(
#   #   list(var_info = data.frame(lake = LAKE, var_info_df),
#   #        sp_richness = data.frame(lake = LAKE, species_richness = nrow(count_dat_example)))
#   # )
#   return(
#     data.frame(lake = LAKE, species_richness = nrow(count_dat_example), var_info_df)
#     )
#   
# }, length_info = zoop_length_total_list0$length_dat_2018, abund_info = count_dat_processed_list$length_dat_2018)
# 
# prac_var_decomp$`Upper Lost`
# 
# prac_var_decomp_df <- do.call(rbind, prac_var_decomp)
# 
# prac_var_decomp_df_removeNA <- prac_var_decomp_df[apply(prac_var_decomp_df, MARGIN = 1, function(x) !any(is.na(x)) ), ]
# 
# prac_var_decomp_df_removeNA %>% 
#   mutate(within_sp_variance = 'within_sp_variance') %>% 
#   ggplot(aes(x = within_sp_variance, y = within_var_prop)) + 
#   geom_violin(fill = '#ffd7ba', alpha = 0.5) +
#   geom_jitter(aes(color = species_richness), width = 0.1, size = 5) +
#   scale_color_gradient(low="#a9d6e5", high="#013a63") +
#   ylim(0, 1.05) +
#   geom_hline(yintercept = c(0, 0.5, 1), color = c('black'), linetype = 'dashed') +
#   theme_bw() +
#   ylab("Within-taxon variance") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#   labs(color = 'Species richness')
#  
# #a9d6e5 --> #013a63
# 
# ength_dat_example <- zoop_length_total_list0$length_dat_2018 %>% 
#   filter(Lake == "Lake Below TB") %>% 
#   select(taxa, length_mm) %>% 
#   rename(species = taxa,
#          trait = length_mm)
# 
# total_count <- count_dat_processed_list$length_dat_2018 %>% 
#   filter(Lake == "Upper Indian Basin") %>% 
#   filter(taxa %in% unique(length_dat_example$species)) %>% 
#   pull(total) %>% 
#   sum()
# 
# count_dat_example <- count_dat_processed_list$length_dat_2018 %>% 
#   filter(Lake == "Upper Indian Basin") %>% 
#   filter(taxa %in% unique(length_dat_example$species)) %>% 
#   mutate(rel_abund = total/total_count) %>% 
#   select(taxa, rel_abund) %>% 
#   rename(species = taxa)
# 
# 
# data.frame(lake = 'x', var_decomp(abundance_df = count_dat_example, trait_df = length_dat_example))
# 
# var_decomp
# trait_info_alt



#################################
### CODE NOT CURRENTLY IN USE ###
#################################


#' library(TPD)
#'
#' count_dat_combine_check <- rbind(count_dat_processed_list$length_dat_2018,
#'      count_dat_processed_list$length_dat_2019) %>%
#'   filter(taxa %in% c("cyclopoid", "cyclopoid  f", "cyclopoid f", "cyclopoid immature", "cyclopoid juvenile", "cyclopoid unident"))
#'
#' sort(unique(count_dat_combine_check$taxa))
#'
#' count_dat_combine <- rbind(count_dat_processed_list$length_dat_2018,
#'                            count_dat_processed_list$length_dat_2019) %>%
#'   mutate(new_lake1 = tolower(gsub(" ", "_", trimws(gsub("Lake", "", Lake))) )) %>%
#'   mutate(new_lake2 = recode(new_lake1, #oldname = newname
#'                             "mid_deep_creek" = "mid_deepcreek",
#'                             "frozen__1" = "frozen1",
#'                             "frozen__2" = "frozen2",
#'                             "low_deep_creek" = "lower_deepcreek",
#'                             "upper_black_joe" = "upper_blackjoe",
#'                             "bear_west" = "west_bear",
#'                             "lower_indian_basin" = "lower_indian"),
#'          new_taxa = recode(taxa,
#'                            "b. longirostris" = "B_longirostris",
#'                            "d. mendatoe" = "D_mendotae",
#'                            "d. middenorfiana" = "D_pulex",
#'                            "d. middenorfiana/pulex" = "D_pulex",
#'                            "h. shoshone" = "H_shoshone",
#'                            "l. minutus" = "L_minutus",
#'                            "cyclopoid" = "UnIdCyclopoid",
#'                            "cyclopoid unident" = "UnIdCyclopoid")) %>%
#'   filter(new_taxa %in% c("B_longirostris", "D_mendotae", "D_pulex", "H_shoshone", "L_minutus","UnIdCyclopoid")) %>%
#'   group_by(new_taxa, new_lake2) %>%
#'   summarize(total_count = sum(total), .groups = "drop") %>%
#'   filter(!is.na(total_count)) %>%
#'   rename(lake = new_lake2,
#'          species = new_taxa)
#'
#' egg_len_dat_processed1 <- egg_len_dat_processed %>%
#'   group_by(species) %>%
#'   mutate(species_size_mm = mean(size_mm),
#'          species_egg = mean(egg_count)) %>%
#'   ungroup() %>%
#'   group_by(species, lake) %>%
#'   summarize(species_lake_size_mm = mean(size_mm),
#'             species_size_mm = first(species_size_mm),
#'             species_lake_egg_count = mean(egg_count),
#'             species_egg = first(species_egg),
#'             contains_fish = first(contains_fish),
#'             .groups = 'drop')
#'
#' abund_data_tpd <- combined_data_tpd[, c('lake', 'contains_fish', 'species_fishinfo', 'count_updated')] %>%
#'   group_by(lake, species_fishinfo) %>%
#'   summarize(contains_fish = first(contains_fish),
#'             count_updated = first(count_updated), .groups = 'drop') %>%
#'   select(!lake) %>%
#'   group_by(species_fishinfo) %>%
#'   summarize(contains_fish = first(contains_fish),
#'             summed_count = sum(count_updated), .groups = 'drop') %>%
#'   pivot_wider(names_from = 'species_fishinfo',
#'               values_from = 'summed_count',
#'               values_fill = 0) %>%
#'   column_to_rownames(var = 'contains_fish') %>%
#'   as.matrix()
#'
#' combined_data_tpd <- left_join(egg_len_dat_processed,
#'                                count_dat_combine,
#'                                by = c('species', 'lake')) %>%
#'   mutate(count_updated = if_else(is.na(total_count), 1, total_count),
#'          species_fishinfo = paste0(species, "_", contains_fish))
#'
#' abund_data_tpd_full_species <- combined_data_tpd[, c('lake', 'contains_fish', 'species', 'count_updated')] %>%
#'   group_by(lake, species) %>%
#'   summarize(contains_fish = first(contains_fish),
#'             count_updated = first(count_updated), .groups = 'drop') %>%
#'   select(!lake) %>%
#'   group_by(species) %>%
#'   summarize(contains_fish = first(contains_fish),
#'             summed_count = sum(count_updated), .groups = 'drop') %>%
#'   pivot_wider(names_from = 'species',
#'               values_from = 'summed_count',
#'               values_fill = 0) %>%
#'   column_to_rownames(var = 'contains_fish') %>%
#'   as.matrix()
#'
#' combined_data_tpd_full_species <- left_join(egg_len_dat_processed,
#'                                count_dat_combine,
#'                                by = c('species', 'lake')) %>%
#'   mutate(count_updated = if_else(is.na(total_count), 1, total_count))
#'
#'
#'
#' #fish vs. fishless traits
#' observed_tpds <- TPDs(species = combined_data_tpd$species_fishinfo,
#'                  traits = combined_data_tpd[, c('size_mm', 'egg_count')])
#' abund_data_tpd_final <- abund_data_tpd[,match(names(observed_tpds$TPDs), colnames(abund_data_tpd))]
#'
#' observed_dissim <- dissim(x = TPDc(TPDs = observed_tpds, abund_data_tpd_final))
#'
#' observed_dissim_df <- data.frame(dissimilarity = observed_dissim$communities$dissimilarity[1,2],
#'            P_shared = observed_dissim$communities$P_shared[1,2],
#'            P_non_shared = observed_dissim$communities$P_non_shared[1,2])
#'
#' #full_species_traits
#' observed_tpds_full_species <- TPDs(species = combined_data_tpd_full_species$species,
#'                       traits = combined_data_tpd_full_species[, c('size_mm', 'egg_count')])
#' abund_data_tpd_final_full_species <- abund_data_tpd_full_species[,match(names(observed_tpds_full_species$TPDs), colnames(abund_data_tpd_full_species))]
#'
# observed_dissim_full_species <- dissim(x = TPDc(TPDs = observed_tpds_full_species, abund_data_tpd_final_full_species))
#
# observed_dissim_df_full_species <- data.frame(dissimilarity = observed_dissim_full_species$communities$dissimilarity[1,2],
#                                  P_shared = observed_dissim_full_species$communities$P_shared[1,2],
#                                  P_non_shared = observed_dissim_full_species$communities$P_non_shared[1,2])
#
#

# before_time <- Sys.time()
# rand_tpd_list <- lapply(1:100, function(x, trait_data, processed_data) {
# 
#     combined_data_tpd_rand <- trait_data %>%
#       group_by(species) %>%
#       #mutate(lake = sample(lake)) %>%
#       mutate(contains_fish = sample(contains_fish)) %>%
#       ungroup() %>%
#       mutate(species_fishinfo = paste0(species, "_", contains_fish))
# 
#   rand_tpd_sp <- TPDs(species = combined_data_tpd_rand$species_fishinfo,
#                       traits = combined_data_tpd_rand[, c('size_mm', 'egg_count')])
# 
#   abund_data_tpd_final <- processed_data[,match(names(rand_tpd_sp$TPDs), colnames(processed_data))]
#   dissim_calc <- dissim(x = TPDc(TPDs = rand_tpd_sp, abund_data_tpd_final))
# 
#   return(
#     data.frame(dissimilarity = dissim_calc$communities$dissimilarity[1,2],
#                P_shared = dissim_calc$communities$P_shared[1,2],
#                P_non_shared = dissim_calc$communities$P_non_shared[1,2])
#   )
# 
# }, trait_data = combined_data_tpd, processed_data = abund_data_tpd)
# 
# rand_tpd_list_collapse <- do.call(rbind, rand_tpd_list)
# 
# after_time <- Sys.time()
# after_time - before_time
# 
# rand_tpd_list_collapse %>%
#   ggplot() +
#   geom_histogram(aes(x = dissimilarity), bins = 20, fill = 'gray') +
#   geom_vline(aes(xintercept = observed_dissim_df$dissimilarity), color = 'black') +
#   theme_bw()
# 
# 
# 
# 
# 
# tpd_s1 <- TPDs(species = rand_tpd_list[[1]]$species_fishinfo,
#      traits = rand_tpd_list[[1]][, c('size_mm', 'egg_count')])
# abund_data_tpd_final1 <- abund_data_tpd[,match(names(tpd_s1$TPDs), colnames(abund_data_tpd))]
# dissim(x = TPDc(TPDs = tpd_s1, abund_data_tpd_final1))
# 
# abund_data_tpd
# tpd_s1
# 
# do.call(rbind, rand_tpd_list)
# 
# 
# 
# 
# 
# 
# 
# traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
# sp_iris <- iris$Species
# TPDs_iris <- TPDs(species = sp_iris, traits_iris)
# 
# names(TPDs_iris$TPDs)
# 
# 
# matrix(c(c(0.5, 0.4, 0), #I. virginica absent
#                                  c(0.0, 0.9,  0.1 ), #I. versic. dominates; setosa absent
#                                  c(0.0, 0.1,  0.9 )), #I. virg. dominates; setosa absent
#                                ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
#                                                                        unique(iris$Species)))
# 
# 
# table(combined_data_tpd$contains_fish,
#       combined_data_tpd$species)
# 
# 
# prac_tpd <- TPDs(species = combined_data_tpd[combined_data_tpd$contains_fish == 'yes',]$species,
#      traits = combined_data_tpd[combined_data_tpd$contains_fish == 'yes', c('size_mm', 'egg_count')])
# 
# combined_data_tpd$contains_fish
# 
# abund_data_tpd <- combined_data_tpd[, c('lake', 'contains_fish', 'species_fishinfo', 'count_updated')] %>%
#   group_by(lake, species_fishinfo) %>%
#   summarize(contains_fish = first(contains_fish),
#             count_updated = first(count_updated), .groups = 'drop') %>%
#   select(!lake) %>%
#   group_by(species_fishinfo) %>%
#   summarize(contains_fish = first(contains_fish),
#             summed_count = sum(count_updated), .groups = 'drop') %>%
#   pivot_wider(names_from = 'species_fishinfo',
#               values_from = 'summed_count',
#               values_fill = 0) %>%
#   column_to_rownames(var = 'contains_fish') %>%
#   as.matrix()
# 
# abund_data_tpd_final <- abund_data_tpd[,match(names(prac_tpd$TPDs), colnames(abund_data_tpd))]
# 
# 
# 
# 
# dissim_calc <- dissim(x = TPDc(TPDs = prac_tpd, abund_data_tpd_final))
# dissim_calc$communities$dissimilarity[1,2]
# 
# data.frame(group = c('a', 'a', 'a', 'b', 'b', 'b', 'c', 'c'),
#            val = 1:8) %>%
#   group_by(group) %>%
#   mutate(val = sample(val))
# 
# 
# abund_data_tpd_final
# 
# 
# 
# abund_data_tpd[,]
# abund_data_tpd[,match(names(prac_tpd$TPDs), colnames(abund_data_tpd))]
# combined_data_tpd$
# 
# 
# TPDs(species = combined_data_tpd[combined_data_tpd$contains_fish == 'no',]$species,
#      traits = combined_data_tpd[combined_data_tpd$contains_fish == 'no', c('size_mm', 'egg_count')])
# 
# 
# 
# TPDs(species = egg_len_dat_processed[egg_len_dat_processed$contains_fish == 'no',]$species,
#      traits = egg_len_dat_processed[egg_len_dat_processed$contains_fish == 'no', c('size_mm', 'egg_count')])
# 
# 
# egg_len_dat_processed[, c('size_mm', 'egg_count')]
# 
# table(egg_len_dat_processed$contains_fish,
#       egg_len_dat_processed$species)
# 
# head(iris)
# 
# 
# count_dat_combine %>%
#   filter(lake == 'lower_windy')
# 
# anti_join(egg_len_dat_processed1,
#            count_dat_combine,
#            by = c('lake', 'species')) %>%
#   as.data.frame()
# 
# 
# egg_len_count_combine <- inner_join(egg_len_dat_processed1,
#                             count_dat_combine,
#                             by = c('lake', 'species'))
# 
# egg_len_count <-  egg_len_count_combine %>%
#   group_by(lake) %>%
#   summarize(species_lake_size_cwm = weighted.mean(species_lake_size_mm, w = total_count),
#             species_size_mm_cwm = weighted.mean(species_size_mm, w = total_count),
#             species_lake_egg_count_cwm = weighted.mean(species_lake_egg_count, w = total_count),
#             species_egg_cwm = weighted.mean(species_egg, w = total_count),
#             contains_fish = first(contains_fish),
#             .groups = 'drop')
# 
# cwm_plot <- egg_len_count %>%
#   ggplot() +
#   geom_point(aes(x = species_lake_size_cwm, y = species_lake_egg_count_cwm, color = contains_fish),
#              size = 5.5, alpha = 0.5) +
#   geom_rug(data = egg_len_count %>%
#              group_by(contains_fish) %>%
#              summarize(mean_size = mean(species_lake_size_cwm),
#                        mean_egg = mean(species_lake_egg_count_cwm)),
#            aes(x = mean_size, y = mean_egg, color = contains_fish),
#            size = 4) +
#   xlab('Size (mm)') + ylab('Egg count') +
#   theme_bw()
# ggsave(filename = here('figures', 'com_specific_cwm_plot_3_29_2022.png'),
#        plot = cwm_plot,
#        width = 12*0.8, height = 8*0.8)
# 
# 
# egg_len_count_combine %>%
#   ggplot() +
#   geom_violin(aes(x = total_count, y = species, color = contains_fish)) +
#   geom_point(aes(x = total_count, y = species, color = contains_fish)) +
#   theme_bw()
# 
# 
# 
# sort(unique(count_dat_combine$taxa))
# sort(unique(egg_len_dat_processed$species))


# sort(unique(count_dat_combine$Lake))
# 
# count_dat_combine %>%
#   mutate(new_lake1 = tolower(gsub(" ", "_", trimws(gsub("Lake", "", Lake))) )) %>%
#   pull(new_lake1) %>%
#   unique() %>%
#   sort()
# sort(unique(egg_len_dat_processed$lake))
# 
# unique(count_dat_combine$taxa)
# 
# unique(egg_len_dat_processed$species)
# 
# ggplot(egg_len_dat_processed, aes(egg_count, fill = contains_fish)) +
#   geom_histogram(bins = 15) +
#   theme_bw() +
#   facet_wrap(vars(species), nrow = 2, scales = 'free') +
#   xlab('Number of eggs') +
#   ylab('Count')
# 
# eggcount_plot <- ggplot() +
#   geom_histogram(data = egg_len_dat_processed %>% filter(contains_fish == 'yes' & lake != 'lucia'),
#                  aes(egg_count, fill = contains_fish),
#                  alpha = 0.4, bins = 10) +
#   geom_histogram(data = egg_len_dat_processed %>% filter(contains_fish == 'no' & lake != 'lucia'),
#                  aes(egg_count, fill = contains_fish),
#                  alpha = 0.4, bins = 10) +
#   theme_bw() +
#   scale_fill_manual(name = NULL,
#                     values = c('blue', 'red'),
#                     guide = guide_legend(override.aes = list(alpha = 0.4) ) ) +
#   facet_wrap(vars(species), nrow = 2, scales = 'free') +
#   xlab('Number of eggs') +
#   ylab('Count')
# 
# ggsave(filename = here('figures', 'egg_count_multipanel_3_29_2022.png'),
#        plot = eggcount_plot,
#        #width = 12*0.8, height = 8*0.8)
# 
# # & lake != 'lucia'
# length_at_repro_plot <- ggplot(egg_len_dat_processed) +
#   geom_histogram(data = egg_len_dat_processed %>% filter(contains_fish == 'yes'),
#                  aes(size_mm, fill = contains_fish),
#                  bins = 10, alpha = 0.4) +
#   geom_histogram(data = egg_len_dat_processed %>% filter(contains_fish == 'no'),
#                  aes(size_mm, fill = contains_fish, alpha = contains_fish),
#                  bins = 10, alpha = 0.4) +
#   scale_fill_manual(name = NULL,
#                     values = c('blue', 'red'),
#                     guide = guide_legend(override.aes = list(alpha = 0.4) ) ) +
#   guides(fill = guide_legend(override.aes = list(alpha = 0.4))) +
#   theme_bw() +
#   facet_wrap(vars(species), nrow = 2, scales = 'free') +
#   xlab('Size at reproduction (mm)') +
#   ylab('Count')
# 
# ggsave(filename = here('figures', 'length_at_repro_multipanel_3_29_2022.png'),
#        plot = length_at_repro_plot,
#        width = 12*0.8, height = 8*0.8)
# 
# 
# eggcount_vs_lengatrepro <- ggplot(egg_len_dat_processed) +
#   geom_point(aes(x = size_mm, y = egg_count, color = contains_fish), size = 3, alpha = 0.5) +
#   stat_smooth(aes(x = size_mm, y = egg_count, group = contains_fish, color = contains_fish),
#               method = "lm", se = FALSE, size = 1.25, linetype = 1) +
#   scale_color_manual(name = NULL, values = c('blue', 'red')) +
#   theme_bw() +
#   facet_wrap(vars(species), scales = 'free', nrow = 2) +
#   xlab('Size at reproduction (mm)') +
#   ylab('Number of eggs')
# 
# ggsave(filename = here('figures', 'eggcount_vs_lengatrepro_multipanel_3_29_2022.png'),
#        plot = eggcount_vs_lengatrepro,
#        width = 12*0.8, height = 8*0.8)
# 
# 
# egg_len_dat_processed %>%
#   group_by(lake) %>%
#   summarize(mean_egg = mean(egg_count),
#             mean_size = mean(size_mm),
#             contains_fish = first(contains_fish)) %>%
#   ggplot() +
#   geom_point(aes(x = mean_egg, y = mean_size, color = contains_fish))
# 
# 
# rbind(zoop_length_total_list0$length_dat_2018[!duplicated(zoop_length_total_list0$length_dat_2018$Lake), c('year', 'Lake')],
#       zoop_length_total_list0$length_dat_2019[!duplicated(zoop_length_total_list0$length_dat_2019$Lake), c('year', 'Lake')]) %>%
#   group_by(Lake) %>%
#   summarize(count = n()) %>%
#   as.data.frame()


# zoop_length_total_list0$length_dat_2019[!duplicated(zoop_length_total_list0$length_dat_2019$Lake), c('year', 'Lake')] %>%
#   group_by(Lake) %>%
#   summarize(count = n()) %>%
#   as.data.frame()
# 
# rbind(zoop_length_total_list0$length_dat_2018,
#       zoop_length_total_list0$length_dat_2019)
# 
# 
# library(TPD)
# 
# count_processed_yearcombine <- lapply(count_dat_processed_list, function(x) {
#   x %>%
#     filter(total != 0 & !is.na(total))
# }) %>%
#   bind_rows()
# 
# 
# count_processed_yearcombine %>%
#   filter(Lake == 'Upper Twin') %>%
#   filter(taxa == 'chydorus')
# 
# zoop_length_total_list %>%
#   filter(Lake == 'Upper Twin') %>%
#   group_by(taxa) %>%
#   summarize(count = n())



#"Upper Spider", "Footprint", "Blackrock" --> cause errors with TPD
#Rapid, Blue, Hidden, Sad Shrimp, Deep --> count data doesn't contain all the taxa that the length data has
# zoop_length_total_list <- zoop_length_total_list0 %>%
#   bind_rows() %>%
#   #filter(!(Lake %in% c("Upper Spider", "Footprint", "Blackrock", "Rapid", "Blue", "Hidden", "Sad Shrimp", "Deep", "Crater", 'Upper Twin'))) %>%
#   group_by(Lake, taxa) %>%
#   filter(n() > 1) %>%
#   ungroup() %>%
#   filter(!is.na(length_mm))

#'Zigzag'
# sort(unique(zoop_length_total_list$Lake))
# 
# sort(unique(count_processed_yearcombine$taxa))
# sort(unique(zoop_length_total_list$taxa))
# 
# count_processed_yearcombine_taxa_subset <- lapply(unique(zoop_length_total_list$Lake), function(lake, count, length) {
# 
#   taxa_in_length <- unique(length[length$Lake == lake,]$taxa)
# 
#   message(lake)
# 
#   processed_count <- count %>%
#     filter((Lake == lake) & (taxa %in% taxa_in_length) )
# 
# 
#   if (!identical(sort(taxa_in_length), sort(unique(processed_count[processed_count$Lake == lake,]$taxa)) ))
#     stop('not identical')
# 
#   return(processed_count)
# 
# }, count = count_processed_yearcombine, length = zoop_length_total_list) %>%
#   bind_rows()



# prac_tpd <- lapply(unique(zoop_length_total_list$Lake)[!(unique(zoop_length_total_list$Lake) %in% c("Upper Spider", "Footprint", "Blackrock") )], function(x) {
#   message(x)
#   TPDs(species = zoop_length_total_list %>%
#          #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#          filter(Lake == x) %>%
#          pull(taxa),
#        traits = zoop_length_total_list %>%
#          #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#          filter(Lake == x) %>%
#          pull(length_mm),
#        samples = zoop_length_total_list %>%
#          #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#          filter(Lake == x) %>%
#          pull(Lake)
#   )
# })

# zoop_length_total_list %>%
#   filter(Lake == 'Upper Spider') %>%
#   as.data.frame() %>%
#   group_split(taxa)



# prac_tpd <- TPDs(species = zoop_length_total_list %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(taxa),
#                  traits = zoop_length_total_list %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(length_mm),
#                  samples = zoop_length_total_list %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(Lake)
# )




# length_processed_yearcombine <- lapply(zoop_length_total_list, function(x) {
#   x %>%
#     filter(total != 0 & !is.na(total))
# }) %>%
#   bind_rows()


# table(zoop_length_total_list$length_dat_2018$Lake,
#       zoop_length_total_list$length_dat_2018$taxa)
# 
# zoop_length_total_list$length_dat_2018 %>%
#   filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#   pull(taxa)

# prac_tpd <- TPDs(species = zoop_length_total_list$length_dat_2018 %>%
#                    group_by(Lake, taxa) %>%
#                    filter(n() > 2) %>%
#                    ungroup() %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(taxa),
#                  traits = zoop_length_total_list$length_dat_2018 %>%
#                    group_by(Lake, taxa) %>%
#                    filter(n() > 2) %>%
#                    ungroup() %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(length_mm),
#                  samples = zoop_length_total_list$length_dat_2018 %>%
#                    group_by(Lake, taxa) %>%
#                    filter(n() > 2) %>%
#                    ungroup() %>%
#                    #filter(Lake %in% c("Little Island", "Mountain Sheep")) %>%
#                    pull(Lake)
#      )



# zoop_length_total_list$length_dat_2018 %>%
#   group_by(Lake, taxa) %>%
#   filter(n() > 2) %>%
#   ungroup()
# 
# example_abundances <- matrix(c(c(0.5, 0.3, 0.2,
#                                  0.1, 0.8, 0.1,
#                                  0.5, 0,   0.5)), #I. virg. dominates; setosa absent
#                              ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
#                                                                      unique(iris$Species)))




# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/eval_func_space_internals_phylofunc.R"))
# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/eval_func_space_phylofunc.R"))
# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/parse_eval_phylofunc.R"))
# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/utilities_phylofunc.R"))
# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/visualization_phylofunc.R"))
# source(here("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Phylogenetic_Community_Project/Scripts/main_functions/calc_fd_devel.R"))



# ### load data ###
# #the rest of the data (length, counts) are loaded and processed together further down in the script
# #egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length.xlsx'), sheet = 'Data'))
# #egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_2_25_2022.xlsx'), sheet = 'Data'))
# egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_3_29_2022.xlsx'), sheet = 'Data formatted for R'))
#
#
#
# #######################
# ### PROCESSING DATA ###
# #######################
#
# ### new processing step on 3/29/2022 ###
# # rename to columns to convenient names ###
# egg_len_dat_processed <- egg_len_dat %>%
#   rename("size_mm" = "Size (mm)",
#          "egg_count" = "Number of eggs") %>%
#   rename_with(.fn = tolower, .cols = everything()) %>%
#   mutate(contains_fish = ifelse(fish == 1, 'yes', 'no'))

# ### 3/29/2022: THIS IS NOT LONGER NECESSARY WITH THE DATA REFORMATTING
# ### Length at reproduction and egg count data ###
# #collected by Spencer Cruz
# #Spencer is going to restructure the data when he is collecting it so that it is in "tidy" format,
# #so much of this processing may no longer be required for the egg count and length at reproduction data
# # egg_len_dat$id <- seq_len(nrow(egg_len_dat)) #add unique ID for each row
# # egg_len_dat$L_minutus_mm[egg_len_dat$L_minutus_mm == 's' & !is.na(egg_len_dat$L_minutus_mm)] <- NA
# # egg_len_dat$L_minutus_mm <- as.numeric(egg_len_dat$L_minutus_mm)
# #
# # egg_len_dat_processed <- egg_len_dat %>%
# #   pivot_longer(cols = B_longirostris_mm:UnIdCyclopoid_eggs, names_to = 'species_trait', values_to = 'value', values_drop_na = TRUE) %>%
# #   mutate(species = gsub("_[A-Za-z]*$", '', species_trait),
# #          trait = gsub("^[A-Za-z]*_*[A-Za-z]*_", '', species_trait),
# #          Fish = if_else(fish == 1, 'fish', 'fishless')) %>%
# #   select(!c(species_trait, fish)) %>%
# #   pivot_wider(names_from = trait, values_from = value) %>%
# #   rename(reproduction_length = mm,
# #          egg_count = eggs)
#
#
# ### LENGTH DATA ###
# #collected by Lindsey Boyle
# #length is of a collection of individuals of each taxa in the lake
# #either all individuals were measured if only a few were in the sample or a random subset of individuals
# #were measured if many individuals were collected (too many to be comprehensively measured)
# zoop_length_total_list0 <- list()
#
# for (YEAR in c(2018, 2019)) {
#
#   wr_zoop_init <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
#     message(NAME)
#     len_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 5, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#
#     #remove all extra columns (everything that is found to the right of the first column labelled as ...NUMBER)
#     #len_dat_1 <- len_dat[,!grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))]
#     extra_col_detect <- grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))
#     if (any(extra_col_detect)) {
#       remove_col_indices <- min(which(extra_col_detect)):ncol(len_dat)
#       len_dat_1 <- len_dat[, !(seq_len(ncol(len_dat)) %in% remove_col_indices) ]
#     } else {
#       len_dat_1 <- len_dat
#     }
#
#     len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(grepl(pattern = 'Averages', len_dat_1[,1])):nrow(len_dat_1)),] #remove the average information
#     #len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(len_dat_1[,1] == 'Averages'):nrow(len_dat_1)),]
#
#     len_dat_3 <- len_dat_2[apply(len_dat_2, MARGIN = 1, function(x) !all(is.na(x))),] #remove any rows with all NAs (e.g. extra rows at the bottom of the dataframe)
#     colnames(len_dat_3) <- tolower(colnames(len_dat_3)) #sometimes "1x" is labelled as "1X", which messes things up
#
#     return(len_dat_3)
#   })
#
#
#   wr_zoop_1 <- lapply(wr_zoop_init, function(DATA) {
#     as.data.frame(apply(DATA, 2, as.numeric)) %>%
#       mutate(id = 1:n()) %>% #add unique id for each row (necessary for correct wide pivot)
#       pivot_longer(cols = !id, names_to = 'taxa_trait', values_to = 'value', values_drop_na = TRUE) %>%
#       mutate(taxa = tolower(str_trim(str_remove(string = taxa_trait, pattern = "1x$|2x$|mm$"))),
#              taxa = ifelse(taxa %in% c("polyphemus"), "p. pediculus", taxa),
#              taxa = gsub("\\.$", "", taxa),
#              measurement = paste0('length_', tolower(str_trim(str_extract(string = taxa_trait, pattern = "1x$|2x$|mm$"))))) %>%
#       select(!taxa_trait) %>%
#       pivot_wider(names_from = measurement, values_from = value) %>%
#       select(!id)
#   })
#
#   zoop_length_total_list0[[paste0('length_dat_', YEAR)]] <- as.data.frame(bind_rows(wr_zoop_1, .id = 'Lake'))
#   zoop_length_total_list0[[paste0('length_dat_', YEAR)]]$year <- YEAR
#
# }
#
#
# ### COUNT DATA ###
# #collected by Lindsey Boyle
# #The number of individuals for each taxa. Some of the taxa have more specific count information such as
# #counts by sex, counts of individuals with eggs, etc...
# count_dat_processed_list <- list()
#
# for (YEAR in c(2018, 2019)) {
#
#   count_dat_processed_list_year <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
#     message(NAME)
#     count_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#     count_dat_1 <- count_dat[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat))] #remove any extra columns (named "...NUMBER")
#
#     count_dat_2 <- count_dat_1 %>%
#       pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#       mutate(taxa_counttype = str_replace(taxa_counttype, pattern = 'juvenile|juveline', replace = 'immature'), #switch juvenile and a mispelling of that to immature
#              count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$"))), #add count type (e.g. m, f, total) to new column
#              taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|[Ff]*[ ]*w/[ ]*eggs$|[Ff]*[ ]*w/[ ]*epiphia$|immature$|\\?")))) %>% #add taxa to new column
#       mutate(taxa = str_remove(taxa, '\\.$|\\?'), #remove trailing periods and question marks from taxa name
#              count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#              count_type = recode(count_type, juvenile = 'immature', juveline = 'immature'),
#              count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#       mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa),
#              taxa = if_else(taxa == "d. mendotae", "d. mendatoe", taxa),
#              taxa = if_else(taxa == "fillinia", "filinia", taxa),
#              taxa = if_else(taxa %in% c("euclanis sp", "euchlanis sp"), "euchlanis", taxa),
#              taxa = if_else(taxa == "h. gibberu", "h. gibberum", taxa),
#              taxa = if_else(taxa == "d. middenorfiana", "d. middenorfiana/pulex", taxa),
#              taxa = ifelse(taxa %in% c("polyphemus pediculus"), "p. pediculus", taxa)) %>%
#       select(!taxa_counttype) %>%
#       pivot_wider(names_from = count_type, values_from = count, values_fill = 0) %>%
#       rowwise() %>%
#       mutate(total = ifelse(total == 0, ifelse(m != 0 | f != 0, sum(m, f, na.rm = TRUE), sum(eggs, na.rm = TRUE)), total)) %>%
#       ungroup() %>%
#       as.data.frame()
#
#     #add immature column if doesnt exist (only one year identified immature/juvenile)
#     if (!('immature' %in% colnames(count_dat_2)))
#       count_dat_2 <- count_dat_2 %>% add_column(immature = 0, .after = "epiphia")
#
#     #replace NAs in the count columns with 0 and remove rows where all counts are 0
#     count_dat_3 <- count_dat_2 %>%
#       replace_na(list(total = 0,
#                       m = 0,
#                       f = 0,
#                       eggs = 0,
#                       immature = 0)
#       ) %>%
#       filter(if_any(c(total, m, f, eggs, immature), ~ . != 0))
#
#     return(count_dat_3)
#   })
#
#   #unique(sapply(count_dat_2018_processed_list, ncol))
#
#   count_dat_processed_list[[paste0('length_dat_', YEAR)]] <- bind_rows(count_dat_processed_list_year, .id = 'Lake')
#   count_dat_processed_list[[paste0('length_dat_', YEAR)]]$year <- YEAR
# }



#fish_info <- egg_len_dat[!duplicated(egg_len_dat$lake), c('lake', 'fish', 'region', 'side')]

# wr_zoop_18_init <- lapply(setNames(nm = excel_sheets(here('data', 'WR2018ZoopsID.xlsx'))), function(NAME) {
#   message(NAME)
#   len_dat_2018 <- as.data.frame(read_excel(here('data', 'WR2018ZoopsID.xlsx'), sheet = NAME, skip = 5, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#
#   #remove all extra columns (everything that is found to the right of the first column labelled as ...NUMBER)
#   #len_dat_2018_1 <- len_dat_2018[,!grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat_2018))]
#   extra_col_detect <- grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat_2018))
#   if (any(extra_col_detect)) {
#     remove_col_indices <- min(which(extra_col_detect)):ncol(len_dat_2018)
#     len_dat_2018_1 <- len_dat_2018[, !(seq_len(ncol(len_dat_2018)) %in% remove_col_indices) ]
#   } else {
#     len_dat_2018_1 <- len_dat_2018
#   }
#
#   len_dat_2018_2 <- len_dat_2018_1[!(seq_len(nrow(len_dat_2018_1)) %in% which(grepl(pattern = 'Averages', len_dat_2018_1[,1])):nrow(len_dat_2018_1)),] #remove the average information
#   #len_dat_2018_2 <- len_dat_2018_1[!(seq_len(nrow(len_dat_2018_1)) %in% which(len_dat_2018_1[,1] == 'Averages'):nrow(len_dat_2018_1)),]
#
#   len_dat_2018_3 <- len_dat_2018_2[apply(len_dat_2018_2, MARGIN = 1, function(x) !all(is.na(x))),] #remove any rows with all NAs (e.g. extra rows at the bottom of the dataframe)
#   colnames(len_dat_2018_3) <- tolower(colnames(len_dat_2018_3)) #sometimes "1x" is labelled as "1X", which messes things up
#
#   return(len_dat_2018_3)
# })
#
#
# wr_zoop_18_1 <- lapply(wr_zoop_18_init, function(DATA) {
#   as.data.frame(apply(DATA, 2, as.numeric)) %>%
#     mutate(id = 1:n()) %>%
#     pivot_longer(cols = !id, names_to = 'taxa_trait', values_to = 'value', values_drop_na = TRUE) %>%
#     mutate(taxa = tolower(str_trim(str_remove(string = taxa_trait, pattern = "1x$|2x$|mm$"))),
#            measurement = paste0('length_', tolower(str_trim(str_extract(string = taxa_trait, pattern = "1x$|2x$|mm$"))))) %>%
#     select(!taxa_trait) %>%
#     pivot_wider(names_from = measurement, values_from = value) %>%
#     select(!id)
# })
#
# wr_zoop_18_collapse <- as.data.frame(bind_rows(wr_zoop_18_1, .id = 'Lake'))
# wr_zoop_18_collapse$year <- 2018
#
#
#
# wr_zoop_19_init <- lapply(setNames(nm = excel_sheets(here('data', 'WR2019ZoopsID.xlsx'))), function(NAME) {
#   message(NAME)
#   len_dat_2019 <- as.data.frame(read_excel(here('data', 'WR2019ZoopsID.xlsx'), sheet = NAME, skip = 5, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#
#   #remove all extra columns (everything that is found to the right of the first column labelled as ...NUMBER)
#   #len_dat_2019_1 <- len_dat_2019[,!grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat_2019))]
#   extra_col_detect <- grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat_2019))
#   if (any(extra_col_detect)) {
#     remove_col_indices <- min(which(extra_col_detect)):ncol(len_dat_2019)
#     len_dat_2019_1 <- len_dat_2019[, !(seq_len(ncol(len_dat_2019)) %in% remove_col_indices) ]
#   } else {
#     len_dat_2019_1 <- len_dat_2019
#   }
#
#   len_dat_2019_2 <- len_dat_2019_1[!(seq_len(nrow(len_dat_2019_1)) %in% which(grepl(pattern = 'Averages', len_dat_2019_1[,1])):nrow(len_dat_2019_1)),] #remove the average information
#   #len_dat_2019_2 <- len_dat_2019_1[!(seq_len(nrow(len_dat_2019_1)) %in% which(len_dat_2019_1[,1] == 'Averages'):nrow(len_dat_2019_1)),]
#
#   len_dat_2019_3 <- len_dat_2019_2[apply(len_dat_2019_2, MARGIN = 1, function(x) !all(is.na(x))),] #remove any rows with all NAs (e.g. extra rows at the bottom of the dataframe)
#   colnames(len_dat_2019_3) <- tolower(colnames(len_dat_2019_3)) #sometimes "1x" is labelled as "1X", which messes things up
#
#   return(len_dat_2019_3)
# })
#
#
# wr_zoop_19_1 <- lapply(wr_zoop_19_init, function(DATA) {
#   as.data.frame(apply(DATA, 2, as.numeric)) %>%
#     mutate(id = 1:n()) %>%
#     pivot_longer(cols = !id, names_to = 'taxa_trait', values_to = 'value', values_drop_na = TRUE) %>%
#     mutate(taxa = tolower(str_trim(str_remove(string = taxa_trait, pattern = "1x$|2x$|mm$"))),
#            measurement = paste0('length_', tolower(str_trim(str_extract(string = taxa_trait, pattern = "1x$|2x$|mm$"))))) %>%
#     select(!taxa_trait) %>%
#     pivot_wider(names_from = measurement, values_from = value) %>%
#     select(!id)
# })
#
# wr_zoop_19_collapse <- as.data.frame(bind_rows(wr_zoop_19_1, .id = 'Lake'))



# count_dat_2018_processed_list <- lapply(setNames(nm = excel_sheets(here('data', 'WR2018ZoopsID.xlsx'))), function(NAME) {
#   message(NAME)
#   count_dat_2018 <- as.data.frame(read_excel(here('data', 'WR2018ZoopsID.xlsx'), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#   count_dat_2018_1 <- count_dat_2018[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat_2018))] #remove any extra columns (named "...NUMBER")
#
#   count_dat_2018_2 <- count_dat_2018_1 %>%
#     pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#     mutate(count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = " [Tt]otal$|M$|F$|w/[ ]*eggs$|w/[ ]*epiphia$"))), #add count type (e.g. m, f, total) to new column
#            taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = " [Tt]otal$|M$|F$|w/[ ]*eggs$|w/[ ]*epiphia$")))) %>% #add taxa to new column
#     mutate(taxa = str_remove(taxa, '\\.$'), #remove trailing periods from taxa name
#            count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#            count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#     mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa)) %>% 
#     select(!taxa_counttype) %>% 
#     pivot_wider(names_from = count_type, values_from = count) %>%
#     as.data.frame()
#   
#   return(count_dat_2018_2)
# })
# 
# #unique(sapply(count_dat_2018_processed_list, ncol))
# 
# count_dat_2018_processed_df <- bind_rows(count_dat_2018_processed_list, .id = 'Lake')
# 
# 
# count_dat_2019_processed_list <- lapply(setNames(nm = excel_sheets(here('data', 'WR2019ZoopsID.xlsx'))), function(NAME) {
#   message(NAME)
#   count_dat_2019 <- as.data.frame(read_excel(here('data', 'WR2019ZoopsID.xlsx'), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#   count_dat_2019_1 <- count_dat_2019[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat_2019))] #remove any extra columns (named "...NUMBER")
#   
#   count_dat_2019_2 <- count_dat_2019_1 %>% 
#     pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#     mutate(count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = " [Tt]otal$|M$|F$|w/[ ]*eggs$|w/[ ]*epiphia$"))), #add count type (e.g. m, f, total) to new column
#            taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = " [Tt]otal$|M$|F$|w/[ ]*eggs$|w/[ ]*epiphia$")))) %>% #add taxa to new column
#     mutate(taxa = str_remove(taxa, '\\.$'), #remove trailing periods from taxa name
#            count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#            count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#     mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa)) %>% 
#     select(!taxa_counttype) %>% 
#     pivot_wider(names_from = count_type, values_from = count) %>%
#     as.data.frame()
#   
#   return(count_dat_2019_2)
# })




# 
# ##########################
# ### SCRIPT PREPARATION ###
# ##########################
# 
# ### load libraries ###
# library(here)
# library(tidyverse)
# library(readxl)
# 
# ### load data ###
# ##the rest of the data (length, counts) are loaded and processed together further down in the script
# ##egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length.xlsx'), sheet = 'Data'))
# ##egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_2_25_2022.xlsx'), sheet = 'Data'))
# 
# # egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_3_29_2022.xlsx'), sheet = 'Data formatted for R'))
# # 
# # 
# # 
# # egg_len_dat_processed <- egg_len_dat %>% 
# #   rename("size_mm" = "Size (mm)",
# #          "egg_count" = "Number of eggs") %>% 
# #   rename_with(.fn = tolower, .cols = everything()) %>% 
# #   mutate(contains_fish = ifelse(fish == 1, 'yes', 'no'))
# 
# 
# ### LENGTH DATA ###
# #collected by Lindsey Boyle
# #length is of a collection of individuals of each taxa in the lake
# #either all individuals were measured if only a few were in the sample or a random subset of individuals
# #were measured if many individuals were collected (too many to be comprehensively measured)
# zoop_length_total_list0 <- list()
# 
# for (YEAR in c(2018, 2019)) {
#   
#   wr_zoop_init <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
#     message(NAME)
#     len_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 5, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#     
#     #remove all extra columns (everything that is found to the right of the first column labelled as ...NUMBER)
#     #len_dat_1 <- len_dat[,!grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))]
#     extra_col_detect <- grepl(pattern = "^\\.{3}[0-9]*", colnames(len_dat))
#     if (any(extra_col_detect)) {
#       remove_col_indices <- min(which(extra_col_detect)):ncol(len_dat)
#       len_dat_1 <- len_dat[, !(seq_len(ncol(len_dat)) %in% remove_col_indices) ]
#     } else {
#       len_dat_1 <- len_dat
#     }
#     
#     len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(grepl(pattern = 'Averages', len_dat_1[,1])):nrow(len_dat_1)),] #remove the average information
#     #len_dat_2 <- len_dat_1[!(seq_len(nrow(len_dat_1)) %in% which(len_dat_1[,1] == 'Averages'):nrow(len_dat_1)),]
#     
#     len_dat_3 <- len_dat_2[apply(len_dat_2, MARGIN = 1, function(x) !all(is.na(x))),] #remove any rows with all NAs (e.g. extra rows at the bottom of the dataframe)
#     colnames(len_dat_3) <- tolower(colnames(len_dat_3)) #sometimes "1x" is labelled as "1X", which messes things up
#     
#     return(len_dat_3)
#   })
#   
#   
#   wr_zoop_1 <- lapply(wr_zoop_init, function(DATA) {
#     as.data.frame(apply(DATA, 2, as.numeric)) %>%
#       mutate(id = 1:n()) %>% #add unique id for each row (necessary for correct wide pivot)
#       pivot_longer(cols = !id, names_to = 'taxa_trait', values_to = 'value', values_drop_na = TRUE) %>% 
#       mutate(taxa = tolower(str_trim(str_remove(string = taxa_trait, pattern = "1x$|2x$|mm$"))),
#              taxa = ifelse(taxa %in% c("polyphemus"), "p. pediculus", taxa),
#              taxa = gsub("\\.$", "", taxa),
#              measurement = paste0('length_', tolower(str_trim(str_extract(string = taxa_trait, pattern = "1x$|2x$|mm$"))))) %>% 
#       select(!taxa_trait) %>% 
#       pivot_wider(names_from = measurement, values_from = value) %>% 
#       select(!id)
#   })
#   
#   zoop_length_total_list0[[paste0('length_dat_', YEAR)]] <- as.data.frame(bind_rows(wr_zoop_1, .id = 'Lake'))
#   zoop_length_total_list0[[paste0('length_dat_', YEAR)]]$year <- YEAR
#   
# }
# 
# 
# ### COUNT DATA ###
# #collected by Lindsey Boyle
# #The number of individuals for each taxa. Some of the taxa have more specific count information such as
# #counts by sex, counts of individuals with eggs, etc...
# count_dat_processed_list <- list()
# 
# for (YEAR in c(2018, 2019)) {
#   
#   count_dat_processed_list_year <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
#     message(NAME)
#     count_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#     count_dat_1 <- count_dat[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat))] #remove any extra columns (named "...NUMBER")
#     
#     count_dat_2 <- count_dat_1 %>% 
#       pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#       mutate(taxa_counttype = str_replace(taxa_counttype, pattern = 'juvenile|juveline', replace = 'immature'), #switch juvenile and a mispelling of that to immature
#              count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$"))), #add count type (e.g. m, f, total) to new column
#              taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|[Ff]*[ ]*w/[ ]*eggs$|[Ff]*[ ]*w/[ ]*epiphia$|immature$|\\?")))) %>% #add taxa to new column
#       mutate(taxa = str_remove(taxa, '\\.$|\\?'), #remove trailing periods and question marks from taxa name
#              count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#              count_type = recode(count_type, juvenile = 'immature', juveline = 'immature'),
#              count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#       mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa),
#              taxa = if_else(taxa == "d. mendotae", "d. mendatoe", taxa),
#              taxa = if_else(taxa == "fillinia", "filinia", taxa),
#              taxa = if_else(taxa %in% c("euclanis sp", "euchlanis sp"), "euchlanis", taxa),
#              taxa = if_else(taxa == "h. gibberu", "h. gibberum", taxa),
#              taxa = if_else(taxa == "d. middenorfiana", "d. middenorfiana/pulex", taxa),
#              taxa = ifelse(taxa %in% c("polyphemus pediculus"), "p. pediculus", taxa)) %>%
#       select(!taxa_counttype) %>% 
#       pivot_wider(names_from = count_type, values_from = count, values_fill = 0) %>%
#       rowwise() %>% 
#       mutate(total = ifelse(total == 0, ifelse(m != 0 | f != 0, sum(m, f, na.rm = TRUE), sum(eggs, na.rm = TRUE)), total)) %>% 
#       ungroup() %>% 
#       as.data.frame()
#     
#     #add immature column if doesnt exist (only one year identified immature/juvenile)
#     if (!('immature' %in% colnames(count_dat_2)))
#       count_dat_2 <- count_dat_2 %>% add_column(immature = 0, .after = "epiphia")
#     
#     #replace NAs in the count columns with 0 and remove rows where all counts are 0
#     count_dat_3 <- count_dat_2 %>% 
#       replace_na(list(total = 0, 
#                       m = 0, 
#                       f = 0, 
#                       eggs = 0, 
#                       immature = 0)
#       ) %>%
#       filter(if_any(c(total, m, f, eggs, immature), ~ . != 0))
#     
#     return(count_dat_3)
#   })
#   
#   #unique(sapply(count_dat_2018_processed_list, ncol))
#   
#   count_dat_processed_list[[paste0('length_dat_', YEAR)]] <- bind_rows(count_dat_processed_list_year, .id = 'Lake')
#   count_dat_processed_list[[paste0('length_dat_', YEAR)]]$year <- YEAR
# }
# 
# 
# 
# count_dat_processed_list$length_dat_2018
# 
# 
# table(count_dat_processed_list$length_dat_2018$Lake),
# count_dat_processed_list$length_dat_2018$taxa)
# 
# 
# 
# data.frame(c1 = letters[1:5],
#            c2 = c(1, 0, 0, 2, 0),
#            c3 = c(1, 0, 3, 9, 0),
#            c4 = c(1, 0, 4, 1, 0)) %>% 
#   filter(if_any(c(c2, c3, c4), ~ . != 0))
# filter(c2 != 0 | c3 != 0 | c4 != 0)
# 
# data.frame(x = c(0, NA, 1, NA),
#            y = c(1, 2, 3, NA),
#            z = 1:4) %>% 
#   replace_na(list(x = 0, y = 0, z = 0))
# 
# 
# add_column(z = -1:1, w = 0)
# 
# #keratella, keratella 2 spine, keratella 1 spine in 2018
# #should d. middenorfiana/pulex and d. middenorfiana be treated separately or be lumped together?
# #difference between immature and juveline?
# 
# #Cyclapoid vs Cyclopoid unident. --> from the same lake
# 
# 
# unique(zoop_length_total_list0$length_dat_2019$taxa)
# unique(count_dat_processed_list$length_dat_2019$taxa)
# unique(count_dat_processed_list$length_dat_2018$taxa)
# 
# unique(zoop_length_total_list0$length_dat_2019$taxa)[]
# 
# count_dat_processed_list$length_dat_2018 %>% 
#   filter(taxa == "keratella")
# 
# unique(count_dat_processed_list$length_dat_2019$taxa)[!unique(count_dat_processed_list$length_dat_2019$taxa) %in% unique(count_dat_processed_list$length_dat_2018$taxa)]
# 
# unique(count_dat_processed_list$length_dat_2018$taxa)[!unique(count_dat_processed_list$length_dat_2018$taxa) %in% unique(count_dat_processed_list$length_dat_2019$taxa)]
# 
# 
# as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = 'Spider', skip = 1, n_max = 1, na = 'NA'))
# 
# 
# data.frame(a = c('abcHELLOdef', 'alksfHELLNOsdf', 'green')) %>% 
#   mutate(a = str_replace(a, pattern = 'HELLO|HELLNO', replace = 'hello'))
# 
# case_when(taxa == "d. mendotae" ~ "d. mendatoe",
#           taxa == "fillinia" ~ "filinia",
#           taxa %in% c("polyphemus pediculus") ~ "p. pediculus",
#           taxa == "h. gibberu" ~ "h. gibberum",
#           !(taxa %in% c("d. mendotae", "fillinia", "polyphemus pediculus", "h. gibberu")) ~ taxa)
# 
# count_dat_processed_list_year <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', 2019, 'ZoopsID.xlsx')))), function(NAME) {
#   message(NAME)
#   count_dat <- as.data.frame(read_excel(here('data', paste0('WR', 2019, 'ZoopsID.xlsx')), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
#   count_dat_1 <- count_dat[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat))] #remove any extra columns (named "...NUMBER")
#   
#   return(count_dat_1 %>% 
#            pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#            mutate(count_type = str_trim(tolower(str_extract(str_trim(taxa_counttype), pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$|\\?"))), #add count type (e.g. m, f, total) to new column
#                   taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$|\\?")))) %>% #add taxa to new column
#            mutate(taxa = str_remove(taxa, '\\.$'), #remove trailing periods from taxa name
#                   count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#                   count_type = recode(count_type, juvenile = 'immature', juveline = 'immature'),
#                   count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#            mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa),
#                   taxa = if_else(taxa == "d. mendotae", "d. mendatoe", taxa),
#                   taxa = if_else(taxa == "fillinia", "filinia", taxa),
#                   taxa = if_else(taxa == "euclanis sp", "euchlanis sp", taxa),
#                   taxa = if_else(taxa == "h. gibberu", "h. gibberum", taxa),
#                   taxa = ifelse(taxa %in% c("polyphemus pediculus"), "p. pediculus", taxa)))
#   count_dat_2 <- count_dat_1 %>% 
#     pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
#     mutate(count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$"))), #add count type (e.g. m, f, total) to new column
#            taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$")))) %>% #add taxa to new column
#     mutate(taxa = str_remove(taxa, '\\.$'), #remove trailing periods from taxa name
#            count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
#            count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
#     mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa),
#            taxa = if_else(taxa == "d. mendotae", "d. mendatoe", taxa),
#            taxa = if_else(taxa == "fillinia", "filinia", taxa),
#            taxa = ifelse(taxa %in% c("polyphemus pediculus"), "p. pediculus", taxa)) %>%
#     select(!taxa_counttype) %>% 
#     pivot_wider(names_from = count_type, values_from = count) %>%
#     rowwise() %>% 
#     mutate(total = ifelse(total == 0, ifelse(m != 0 | f != 0, sum(m, f, na.rm = TRUE), sum(eggs, na.rm = TRUE)), total)) %>% 
#     ungroup() %>% 
#     as.data.frame()
#   
#   return(count_dat_2)
# })
# 
# 
# data.frame(taxa_counttype = c("cyclopoid  f", "cyclopoid f", "syncheata ?")) %>%
#   mutate(count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$|\\?"))), #add count type (e.g. m, f, total) to new column
#          taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$|juvenile$|juveline$|\\?"))))
# 
# lapply(count_dat_processed_list_year , function(x) unique(x$taxa))
# 
# as.data.frame(count_dat_processed_list_year$Bluebell)
# 
# 
# 
# 
# immature|juvenile|juveline|
#   
#   
#   #fillinia --> filinia
#   
#   "h. shoshone juveline" --> "h. shoshone juvenile"
# "syncheata ?" --> "syncheata"
# "euchlanis sp" --> "euchlanis"
# "euclanis sp" --> "euchlanis"
# 
# "cyclopoid juvenile" --> "cyclopoid immature"

# count_dat_processed_list_name_update <- lapply(count_dat_processed_list, function(x) {
#   x %>% 
#     mutate(Lake = recode(Lake, "north_blue" = "n of blue"))
# })


### PROCESSING MOVED TO PROCESSING SCRIPTS
# #taxa included in initial analysis
# focal_taxa <- c('b. longirostris',
#                 'd. mendatoe',
#                 'd. middenorfiana/pulex',
#                 "h. shoshone",
#                 'l. minutus',
#                 'chydorus',
#                 'h. gibberum',
#                 'p. pediculus')
# 
# ### PROCESSING TRAIT DATA ###
# trait_info <- trait_info_init %>% 
#   select(Species, `Reproductive Mode`, `Body Shape`, `Feeding Type`) %>% 
#   rename(species = Species,
#          reproductive_mode = `Reproductive Mode`,
#          body_shape = `Body Shape`,
#          feeding_type = `Feeding Type`)
# 
# length_dat_processed <- zoop_length_total_list0 %>%
#   bind_rows() %>% 
#   group_by(taxa) %>% 
#   summarise(mean_length = mean(length_mm, na.rm = TRUE), .groups = 'drop')
# 
# trait_info <- trait_info %>% 
#   na.omit() %>% 
#   mutate(taxa = case_when(species == 'B_longirostris' ~ 'b. longirostris',
#                           species == 'D_mendotae' ~ 'd. mendatoe',
#                           species == 'D_pulex' ~ 'd. middenorfiana/pulex',
#                           species == 'H_shoshone' ~ "h. shoshone",
#                           species == 'L_minutus' ~ 'l. minutus',
#                           species == 'Chydorus' ~ 'chydorus',
#                           species == 'H. gibberum' ~ 'h. gibberum',
#                           species == 'P. pediculous' ~ 'p. pediculus')) %>%
#   select(!species) %>% 
#   left_join(., length_dat_processed, by = 'taxa') %>% 
#   column_to_rownames(var = "taxa")

