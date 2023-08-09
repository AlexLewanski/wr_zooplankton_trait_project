#######################################################
#######################################################
### ANLAYSIS OF ZOOPLANKTON DATA --- BETA DIVERSITY ###
#######################################################
#######################################################

#https://github.com/ShanKothari/DecomposingFD
#https://www.biorxiv.org/content/10.1101/530782v1
#https://stats.stackexchange.com/questions/92542/how-to-perform-a-bootstrap-test-to-compare-the-means-of-two-samples


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
#library(picante) #matrix randomization
library(ade4)
library(vegan)

#load custom functions
source(here("scripts", "custom_functions_zooplankton_project.R"))


### LOADING DATA ###
#(processed) count and length data
count_dat_processed_list <- readRDS(file = here('data', 'processed_data', 'count_dat_processed_list.rds'))
zoop_length_total_list0 <- readRDS(file = here('data', 'processed_data', 'zoop_length_total_list0.rds'))

#fish presence/absence info
lindsey_master_data <- as.data.frame(read_excel(here('data', 'WR_data_ALL_Boyle.xlsx')))

#trait data (not including length data)
trait_info_init <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_5_9_2022.xlsx'), sheet = 'Inter Traits'))

all_taxa_traits <- readRDS(here('data', 'processed_data', 'all_taxa_trait_info_processed.rds'))



#########################################################
### PROCESSING OF DATA AND INITIAL STEPS FOR FD CALCS ###
#########################################################

###  FISH PRESENCE/ABSENCE INFO FOR EACH LAKE ###
fish_info <- lindsey_master_data %>% 
  select(lake, fish)

processed_lake_info <- fish_info[!duplicated(fish_info$lake),] %>%
  mutate(lake = recode(lake,
                       'lake37' = "lake_37",
                       'blackjoe' = "black_joe",
                       'upper_blackjoe' = "upper_black_joe",
                       "below_tb" = "lake_below_tb"))


all_taxa_traits_processed <- all_taxa_traits %>% 
  filter(!(species %in% c('euchlanis', 'filinia')) ) %>%
  column_to_rownames(var = 'species') 


### STEP 1: CALCULATE DISSIMILARITY MATRIX ###
#de Bello et al.'s (2021) extension of Gower distance used for dissim. mat. construction
zoo_dist <- dist_wrapper(df = all_taxa_traits_processed,
                         method = c('dist', 'gowdis', 'gawdis')[3], arg_list = NULL) 


### STEP 2: DETERMINE OPTIMAL FUNCTIONAL SPACE ###
zoo_msd_calc <- eval_func_space(trait_data = zoo_dist,
                                dim = 10,
                                eval_option = c("corank", "mad")[1],
                                show_messages = TRUE,
                                rand_rep = 99)

#information about the optimal functional space
### *** WITH THE FINAL DATASET, THIS SHOULD NOT BE TRUNCATED AT 2 DIMENSIONS
optimal_func_space_info <- zoo_msd_calc$func_eval$corank[which(zoo_msd_calc$func_eval$corank$dim <= 4),] %>% 
  filter(auc == max(auc))


count_dat_processed_beta <- count_dat_processed_list %>%
  bind_rows() %>% 
  mutate(Lake = recode(Lake, 
                       "N of Blue" = "North Blue")) %>% 
  filter(Lake != "Lost Lake") %>% #remove lost lake
  group_by(Lake) %>% 
  filter(year == max(year)) %>% 
  ungroup() %>% 
  group_by(Lake, taxa) %>%
  summarize(total_abund = sum(total), .groups = 'drop') %>% 
  filter(taxa %in% unique(rownames(all_taxa_traits_processed) )) %>%
  #group_by(Lake) %>% 
  #filter(n() >= (optimal_func_space_info$dim + 1) ) %>% 
  #ungroup() %>%
  pivot_wider(names_from = Lake, values_from = total_abund, values_fill = 0) %>%
  column_to_rownames(var = 'taxa') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'lake') %>%
  mutate(lake = gsub(" ", "_", str_trim(tolower(lake)))) %>%
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
                       "bap_lake" = 'bap',
                       "cutthroat_(stough)" = 'cutthroat_stough'))

#convert counts to relative abundances
count_dat_processed_beta[,-1] <- t(apply(count_dat_processed_beta[, -1], 1, function(x) x/sum(x)))



####################################################
### BETA DIVERSITY BASED ON SCHEINER ET AL. 2016 ###
####################################################

#all pairwise combination of lakes
pairwise_lake_combos <- as.data.frame(t(combn(count_dat_processed_beta$lake, 2)))
colnames(pairwise_lake_combos) <- c('lake1', 'lake2')

#calculate distances in the final functional space (same as what is used in all other FD calculations)
zoo_dist_processed_func_space <- dist(zoo_msd_calc$supplementary_info$PCoA_list[[optimal_func_space_info$correction]][,1:optimal_func_space_info$dim])


### PAIRWISE BETA DIVERSITY CALCULATIONS ###
beta_progress <- txtProgressBar(min = 0, max = nrow(pairwise_lake_combos), 
                                              char = "*", style = 3)

pairwise_beta_list <- list() #initiate list to store results

for (i in seq_len(nrow(pairwise_lake_combos))) {
  
  lake1 <- pairwise_lake_combos[i, 'lake1'] #first lake
  lake2 <- pairwise_lake_combos[i, 'lake2'] #second lake
  
  #subset count data down to focal lakes (& move lake column to rowname)
  count_dat_processed_beta_subset1 <- count_dat_processed_beta %>% 
    filter(lake %in% c(lake1, lake2) ) %>% 
    column_to_rownames(var = 'lake')
  
  #change order of taxa to match the order of taxa in the distance matrix
  count_dat_processed_beta_subset2 <- count_dat_processed_beta_subset1[, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(count_dat_processed_beta_subset1) )]
  
  #calculate beta diversity for the two lakes
  beta_calc <- FTD.beta(tdmat = zoo_dist_processed_func_space,
                        spmat = count_dat_processed_beta_subset2,
                        abund = TRUE, q = 1)

  #add the beta diversity info (+ info on compared lakes) as a dataframe to the results list
  pairwise_beta_list[[paste(lake1, lake2, sep = "_")]] <- data.frame(lake1 = lake1,
                                                                     lake2 = lake2,
                                                                     lake1_fish = processed_lake_info[processed_lake_info$lake == lake1, 'fish'],
                                                                     lake2_fish = processed_lake_info[processed_lake_info$lake == lake2, 'fish'],
                                                                     as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight']))
  )
  
  setTxtProgressBar(beta_progress, i) #update progress bar
}


pairwise_beta_df <- bind_rows(pairwise_beta_list, .id = 'lake')

count_dat_processed_beta_comdist <- count_dat_processed_beta %>% 
  column_to_rownames(var = 'lake') %>% 
  t()

count_dat_processed_beta_comdist_reorder <- count_dat_processed_beta_comdist[match(rownames(zoo_msd_calc$supplementary_info$PCoA_list[[optimal_func_space_info$correction]][,1:optimal_func_space_info$dim]),
                                                                                   rownames(count_dat_processed_beta_comdist)),]


zooplankton_weighted_mean_list <- list()
for (i in 1:ncol(count_dat_processed_beta_comdist_reorder)) {
  
  #mean of each trait in each community weighted by the abundance of taxa in the community
  zooplankton_weighted_mean_list[[i]] <- as.data.frame(as.list(apply(zoo_msd_calc$supplementary_info$PCoA_list[[optimal_func_space_info$correction]][,1:optimal_func_space_info$dim], 
                                                         2, function(x) weighted.mean(x, w = count_dat_processed_beta_comdist_reorder[,i]))))
  
  zooplankton_weighted_mean_list[[i]]$community <- colnames(count_dat_processed_beta_comdist_reorder)[i]
}

zooplankton_weighted_mean_df <- do.call(rbind, zooplankton_weighted_mean_list) %>% 
  column_to_rownames(var = 'community') 

zooplankton_weighted_mean_dist <- dist(zooplankton_weighted_mean_df)

### PERFORM PERMANOVA ANALYSES ###
# DO THE FISH AND FISHLESS LAKES CLUSTER IN DIFFERENT PLACES IN FUNCTONAL SPACE?

#load the beta results if they don't exist yet
if (!exists('pairwise_beta_df'))
  pairwise_beta_df <- read.csv(here('results', 'pairwise_beta_fd_scheiner.csv'))

#convert dataframe of distances to a matrix
dist_qdtm_beta <- df2dist_update(as.data.frame(pairwise_beta_df[, c('lake1', 'lake2', "qDTM.beta")]) )
dist_qdtm_beta_reorder <- reorder_dist(dist_mat = dist_qdtm_beta, new_order = processed_lake_info$lake)

beta_dist_list <- list(qdtm = dist_qdtm_beta,
                       weighted_mean = zooplankton_weighted_mean_dist)

multidim_analyses <- lapply(beta_dist_list, function(x, lake_info) {
  
  output_list <- list()
  
  output_list[['reorded_dist']] <- reorder_dist(dist_mat = x, 
                                                new_order = lake_info$lake)
  
  output_list[['permanova']] <- adonis2(output_list[['reorded_dist']] ~ fish, 
                                        data = lake_info %>% mutate(fish = as.character(fish)), 
                                        permutations = 9999)
  
  output_list[['permdisp']] <- betadisper(d = output_list[['reorded_dist']], 
                                          group = lake_info %>% mutate(fish = as.factor(fish)) %>% pull(fish), 
                                           type = c("median","centroid")[2], 
                                          bias.adjust = TRUE)
  
  output_list[['permdisp_sig']] <- TukeyHSD(x = output_list[['permdisp']], 
                                             which = "group", ordered = FALSE,
                                             conf.level = 0.95)
  
  return(output_list)
  
}, lake_info = processed_lake_info)


#how similar are the two beta diversity measures?
mantel_beta <- mantel.rtest(multidim_analyses$qdtm$reorded_dist,
                            multidim_analyses$weighted_mean$reorded_dist)

#permanova results
multidim_analyses$qdtm$permanova
multidim_analyses$weighted_mean$permanova

multidim_analyses$qdtm$permdisp
multidim_analyses$qdtm$permdisp_sig

multidim_analyses$weighted_mean$permdisp
multidim_analyses$weighted_mean$permdisp_sig


### VISUALIZE BETA DISTANCE ###
qdtm_beta_pco <- dudi.pco(quasieuclid(dist_qdtm_beta_reorder), scannf = FALSE, full = TRUE)
weighted_com_dist_pco <- dudi.pco(quasieuclid(zooplankton_weighted_mean_dist), scannf = FALSE, full = TRUE)

beta_div_pcoa <- rbind(qdtm_beta_pco$li[1:4] %>% 
        mutate(beta_type = 'qDTM beta') %>% 
        rownames_to_column(var = 'lake'),
      weighted_com_dist_pco$li %>% 
        mutate(beta_type = 'Weighted mean distance') %>% 
        rownames_to_column(var = 'lake')) %>%
  left_join(., 
            processed_lake_info %>% mutate(fish = as.character(fish)), by = 'lake') %>%
  mutate(fish_info = if_else(fish == '1', 'fish', 'fishless')) %>% 
  ggplot() +
  geom_point(aes(x = A1, y = A2, color = fish_info), size = 5, alpha = 0.7) +
  scale_color_manual(values = c('#0A9396', '#BB3E03')) +
  #xlab(paste0('PC1 (', round((weighted_com_dist_pco$eig/sum(weighted_com_dist_pco$eig))*100, 2)[1], '%)') ) + 
  #ylab(paste0('PC2 (', round((weighted_com_dist_pco$eig/sum(weighted_com_dist_pco$eig))*100, 2)[2], '%)') ) + 
  xlab('PC1') + ylab('PC2') +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 15)) +
  facet_wrap(~beta_type, scales = 'free')


### EXPORTING RESULTS ###
#combine scheiner and weighted means into a single dataframe
beta_output_df <- as.data.frame(as.matrix(zooplankton_weighted_mean_dist)) %>% 
  rownames_to_column(var = 'lake1') %>% 
  pivot_longer(cols = !lake1, 
               names_to = 'lake2', 
               values_to = 'weighted_mean') %>% 
  rowwise() %>% 
  mutate(lake_sort = paste0(sort(c(lake1, lake2)), collapse = '_') ) %>% 
  ungroup() %>% 
  filter(lake1 != lake2) %>% 
  left_join(., pairwise_beta_df %>% 
              rowwise() %>% 
              mutate(lake_sort = paste0(sort(c(lake1, lake2)), collapse = '_')) %>% 
              ungroup()
  ) %>% 
  select(!c(lake_sort, lake))

write.csv(x = beta_output_df,
          file = here('results', 'pairwise_beta_fd_scheiner.csv'), 
          row.names = FALSE)

ggsave(filename = here('figures', 'beta_div_pcoa_fig.png'), 
       plot = beta_div_pcoa,
       width = 13*0.65, height = 6*0.65, bg = 'white')



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# #permanova
# qdtm.div <- adonis2(dist_qdtm_beta_reorder ~ fish, 
#                     data = processed_lake_info %>% mutate(fish = as.character(fish)), 
#                     permutations = 9999)

# qdtm_beta_pco$li %>% 
#   select(1:5) %>% 
#   rownames_to_column(var = 'lake') %>% 
#   left_join(., 
#           processed_lake_info %>% mutate(fish = as.character(fish)), by = 'lake') %>%
#   mutate(fish_info = if_else(fish == '1', 'fish', 'fishless')) %>% 
#   ggplot() +
#   geom_point(aes(x = A1, y = A2, color = fish_info), size = 5, alpha = 0.7) +
#   scale_color_manual(values = c('#0A9396', '#BB3E03')) +
#   xlab(paste0('PC1 (', round((qdtm_beta_pco$eig/sum(qdtm_beta_pco$eig))*100, 2)[1], '%)') ) + 
#   ylab(paste0('PC2 (', round((qdtm_beta_pco$eig/sum(qdtm_beta_pco$eig))*100, 2)[2], '%)') ) + 
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.title = element_text(size = 15))
#   
# ggsave(filename = here('figures', 'beta_pairwise_pcoa.png'), 
#        width = 13*0.65, height = 9*0.65, bg = 'white')



# weighted_com_dist_pco$li %>% 
#   rownames_to_column(var = 'lake') %>% 
#   left_join(., 
#             processed_lake_info %>% mutate(fish = as.character(fish)), by = 'lake') %>%
#   mutate(fish_info = if_else(fish == '1', 'fish', 'fishless')) %>% 
#   ggplot() +
#   geom_point(aes(x = A1, y = A2, color = fish_info), size = 5, alpha = 0.7) +
#   scale_color_manual(values = c('#0A9396', '#BB3E03')) +
#   xlab(paste0('PC1 (', round((weighted_com_dist_pco$eig/sum(weighted_com_dist_pco$eig))*100, 2)[1], '%)') ) + 
#   ylab(paste0('PC2 (', round((weighted_com_dist_pco$eig/sum(weighted_com_dist_pco$eig))*100, 2)[2], '%)') ) + 
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.title = element_text(size = 15))


# processed_lake_info1 <- processed_lake_info %>% 
#   mutate(fish_info = if_else(fish == 1, 'fish', 'fishless')) %>% 
#   select(lake, fish_info)
# 
# fish_fishless_list <- split(processed_lake_info1, 
#                             processed_lake_info1$fish_info)
# 
# 
# fish_fishless_beta_analysis_list <- lapply(fish_fishless_list, function(x, func_space, count_dat, full_count_data) {
#   
#   lake_vec <- x$lake
#   
#   count_dat_subset <- count_dat %>% 
#     filter(lake %in% lake_vec ) %>% 
#     column_to_rownames(var = 'lake')
#   
#   count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(func_space)), colnames(count_dat_subset) )]
#   
#   
#   beta_calc <- FTD.beta(tdmat = func_space,
#                         spmat = count_dat_subset2,
#                         abund = TRUE, q = 1)
#   
#   
#   rando_subsets_list <- replicate(1000, count_dat_processed_beta %>%
#                                     column_to_rownames(var = 'lake') %>% 
#                                     slice_sample(n = length(lake_vec), replace = TRUE),
#                                   simplify = FALSE)
#   
#   message('starting calculations of random lake subsets for the ', x$fish_info[1], ' lakes')
#   beta_rand_progress <- txtProgressBar(min = 0, max = length(rando_subsets_list), 
#                                        char = "*", style = 3)
#   rand_beta_list <- list()
#   
#   for (i in seq_along(rando_subsets_list)) {
#     
#     count_rand_subset <- rando_subsets_list[[i]][, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(rando_subsets_list[[i]]) )]
#     
#     beta_calc_rand_subset <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                                       spmat = count_rand_subset,
#                                       abund = TRUE, q = 1)
#     
#     rand_beta_list[[paste('rand', i, sep = '_')]] <- as.data.frame((beta_calc_rand_subset[names(beta_calc_rand_subset) != 'disp.mat.weight']))
#     setTxtProgressBar(beta_rand_progress, i)
#     
#   }
#   
#   observed_beta <- data.frame(fish_info = x$fish_info[1],
#                               dataset = 'observed',
#                               as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight'])))
#   
#   rand_beta <- bind_rows(rand_beta_list, .id = 'dataset') %>% 
#     mutate(fish_info = x$fish_info[1]) %>% 
#     relocate(fish_info, .before = dataset)
#   
#   return(list(observed = observed_beta,
#               randomized = rand_beta ))
#   
# }, func_space = zoo_dist_processed_func_space, count_dat = count_dat_processed_beta)
# 
# saveRDS(fish_fishless_beta_analysis_list,
#         here('results', 'fish_fishless_beta_comparetorandom_list.rds'))
# 
# 
# fish_fishless_beta_analysis_list <- readRDS(here('results', 'fish_fishless_beta_comparetorandom_list.rds'))
# 
# 
# fish_fishless_beta_analysis_list$fish$randomized %>% 
#   ggplot() +
#   geom_histogram(aes(x = qDTM.beta), fill = 'gray') +
#   geom_vline(xintercept = fish_fishless_beta_analysis_list$fish$observed$qDTM.beta) +
#   theme_bw()
# 
# 
# fish_fishless_beta_analysis_list$fishless$randomized %>% 
#   ggplot() +
#   geom_histogram(aes(x = qDTM.beta), fill = 'gray') +
#   geom_vline(xintercept = fish_fishless_beta_analysis_list$fishless$observed$qDTM.beta) +
#   theme_bw()
# 
# 
# 
# 
# fish_fishless_df <- bind_rows(fish_fishless_list)
# 
# fish_fishless_df <- rbind(fish_fishless_list$fish[sample(1:nrow(fish_fishless_list$fishless)),], fish_fishless_list$fishless)
# 
# 
# observed_beta_list <- list()
# for (i in c('fish', 'fishless')) {
#     
#   lake_vec <- fish_fishless_df[fish_fishless_df$fish_info == i,]$lake
#     
#   count_dat_subset <- count_dat_processed_beta %>% 
#     filter(lake %in% lake_vec) %>% 
#     column_to_rownames(var = 'lake')
#     
#   count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(count_dat_subset) )]
#     
#   beta_calc <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                         spmat = count_dat_subset2,
#                         abund = TRUE, q = 1)
#     
#   observed_beta_list[[i]] <- as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight']))
#     
# }
#   
#   
# observed_beta_dif <- data.frame(M.beta_diff = observed_beta_list$fishless$M.beta - observed_beta_list$fish$M.beta,
#                                         M.beta.prime_dif = observed_beta_list$fishless$M.beta.prime - observed_beta_list$fish$M.beta.prime,
#                                         Ht.beta_dif = observed_beta_list$fishless$Ht.beta - observed_beta_list$fish$Ht.beta,
#                                         qDT.beta_dif = observed_beta_list$fishless$qDT.beta - observed_beta_list$fish$qDT.beta,
#                                         qDTM.beta_dif = observed_beta_list$fishless$qDTM.beta - observed_beta_list$fish$qDTM.beta)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# permuted_fish_info <- lapply(1:100, function(x, dat) {
#   dat$fish_info <- sample(dat$fish_info)
#   return(dat)
# }, dat = fish_fishless_df)
# 
# 
# beta_rand_progress <- txtProgressBar(min = 0, max = length(permuted_fish_info), 
#                                      char = "*", style = 3)
# permuted_beta_list <- list()
# for (x in 1:length(permuted_fish_info)) {
#   beta_list <- list()
#   for (i in c('fish', 'fishless')) {
#     
#     lake_vec <- permuted_fish_info[[x]][permuted_fish_info[[x]]$fish_info == i,]$lake
#     
#     count_dat_subset <- count_dat_processed_beta %>% 
#       filter(lake %in% lake_vec) %>% 
#       column_to_rownames(var = 'lake')
#     
#     count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(count_dat_subset) )]
#     
#     beta_calc <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                           spmat = count_dat_subset2,
#                           abund = TRUE, q = 1)
#     
#     beta_list[[i]] <- as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight']))
#     
#   }
#   
#   
#   permuted_beta_list[[x]] <- data.frame(M.beta_diff = beta_list$fishless$M.beta - beta_list$fish$M.beta,
#              M.beta.prime_dif = beta_list$fishless$M.beta.prime - beta_list$fish$M.beta.prime,
#              Ht.beta_dif = beta_list$fishless$Ht.beta - beta_list$fish$Ht.beta,
#              qDT.beta_dif = beta_list$fishless$qDT.beta - beta_list$fish$qDT.beta,
#              qDTM.beta_dif = beta_list$fishless$qDTM.beta - beta_list$fish$qDTM.beta)
#   
#   setTxtProgressBar(beta_rand_progress, x)
# }
# 
# fish_fishless_beta_analysis_list
# 
# bind_rows(permuted_beta_list) %>% 
#   ggplot() +
#   geom_histogram(aes(x = qDTM.beta_dif), bins = 15, fill = 'gray') +
#   geom_vline(xintercept = observed_beta_dif$qDTM.beta_dif, color = 'red') +
#   theme_bw()
# 
# permuted_beta_list
# 
# #fishless - fish
# #decreasing difference: the two lakes types have more similar beta diversity
# 
# count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(count_dat_subset) )]
# 
# beta_calc <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                       spmat = count_dat_subset2,
#                       abund = TRUE, q = 1)
# 
# 
# data.frame(fish_info = fish_fishless_list$fish$fish_info[1],
#            dataset = 'observed',
#            as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight'])))
# 
# 
# within_treatment_beta_analysis_list <- lapply(fish_fishless_list, function(x, func_space, count_dat, full_count_data) {
#   
#   lake_vec <- x$lake
#   
#   count_dat_subset <- count_dat %>% 
#     filter(lake %in% lake_vec ) %>% 
#     column_to_rownames(var = 'lake')
# 
#     
#   count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(func_space)), colnames(count_dat_subset) )]
#   
#   beta_calc <- FTD.beta(tdmat = func_space,
#                         spmat = count_dat_subset2,
#                         abund = TRUE, q = 1)
#   
#   
#   rando_subsets_list <- replicate(100, count_dat_subset %>%
#                                     #column_to_rownames(var = 'lake') %>% 
#                                     slice_sample(n = nrow(count_dat_subset), replace = TRUE),
#                                   simplify = FALSE)
#   
#   message('starting calculations of random lake subsets for the ', x$fish_info[1], ' lakes')
#   beta_rand_progress <- txtProgressBar(min = 0, max = length(rando_subsets_list), 
#                                        char = "*", style = 3)
#   rand_beta_list <- list()
#   
#   for (i in seq_along(rando_subsets_list)) {
#     
#     count_rand_subset <- rando_subsets_list[[i]][, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(rando_subsets_list[[i]]) )]
#     
#     beta_calc_rand_subset <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                                       spmat = count_rand_subset,
#                                       abund = TRUE, q = 1)
#     
#     rand_beta_list[[paste('rand', i, sep = '_')]] <- as.data.frame((beta_calc_rand_subset[names(beta_calc_rand_subset) != 'disp.mat.weight']))
#     setTxtProgressBar(beta_rand_progress, i)
#     
#   }
#   
#   observed_beta <- data.frame(fish_info = x$fish_info[1],
#                               dataset = 'observed',
#                               as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight'])))
#   
#   rand_beta <- bind_rows(rand_beta_list, .id = 'dataset') %>% 
#     mutate(fish_info = x$fish_info[1]) %>% 
#     relocate(fish_info, .before = dataset)
#   
#   return(list(observed = observed_beta,
#               randomized = rand_beta ))
#   
# }, func_space = zoo_dist_processed_func_space, count_dat = count_dat_processed_beta)
# 
# saveRDS(within_treatment_beta_analysis_list,
#         here('results', 'within_treatment_beta_analysis_list.rds'))
# 
# within_treatment_beta_analysis_list <- readRDS(here('results', 'within_treatment_beta_analysis_list.rds'))
# 
# 
# 
# 
# rbind(within_treatment_beta_analysis_list$fish$randomized,
#       within_treatment_beta_analysis_list$fishless$randomized) %>% 
#   ggplot() +
#   geom_histogram(aes(x = qDTM.beta, fill = fish_info), position="identity", alpha = 0.5) +
#   theme_bw()
# 
# 
# 
# 
# 
# 
# within_treatment_beta_analysis_largersubset_list <- lapply(fish_fishless_list, function(x, min_sample_size, func_space, count_dat, full_count_data) {
#   
#   lake_vec <- x$lake
#   
#   count_dat_subset <- count_dat %>% 
#     filter(lake %in% lake_vec ) %>% 
#     column_to_rownames(var = 'lake')
#   
#   count_dat_subset2 <- count_dat_subset[, match(rownames(as.matrix(func_space)), colnames(count_dat_subset) )]
#   
#   
#   beta_calc <- FTD.beta(tdmat = func_space,
#                         spmat = count_dat_subset2,
#                         abund = TRUE, q = 1)
#   
#   
#   rando_subsets_list <- replicate(1000, count_dat_subset %>%
#                                     #column_to_rownames(var = 'lake') %>% 
#                                     slice_sample(n = min_sample_size, replace = TRUE),
#                                   simplify = FALSE)
#   
#   message('starting calculations of random lake subsets for the ', x$fish_info[1], ' lakes')
#   beta_rand_progress <- txtProgressBar(min = 0, max = length(rando_subsets_list), 
#                                        char = "*", style = 3)
#   rand_beta_list <- list()
#   
#   for (i in seq_along(rando_subsets_list)) {
#     
#     count_rand_subset <- rando_subsets_list[[i]][, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(rando_subsets_list[[i]]) )]
#     
#     beta_calc_rand_subset <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                                       spmat = count_rand_subset,
#                                       abund = TRUE, q = 1)
#     
#     rand_beta_list[[paste('rand', i, sep = '_')]] <- as.data.frame((beta_calc_rand_subset[names(beta_calc_rand_subset) != 'disp.mat.weight']))
#     setTxtProgressBar(beta_rand_progress, i)
#     
#   }
#   
#   observed_beta <- data.frame(fish_info = x$fish_info[1],
#                               dataset = 'observed',
#                               as.data.frame((beta_calc[names(beta_calc) != 'disp.mat.weight'])))
#   
#   rand_beta <- bind_rows(rand_beta_list, .id = 'dataset') %>% 
#     mutate(fish_info = x$fish_info[1]) %>% 
#     relocate(fish_info, .before = dataset)
#   
#   return(list(observed = observed_beta,
#               randomized = rand_beta ))
#   
# }, func_space = zoo_dist_processed_func_space, count_dat = count_dat_processed_beta, min_sample_size = min(table(processed_lake_info1$fish_info)))
# 
# 
# saveRDS(within_treatment_beta_analysis_largersubset_list,
#         here('results', 'within_treatment_beta_analysis_largersubset_list.rds'))
# 
# within_treatment_beta_analysis_largersubset_list <- readRDS(here('results', 'within_treatment_beta_analysis_largersubset_list.rds'))
# 
#       
# rbind(within_treatment_beta_analysis_largersubset_list$fish$randomized,
#       within_treatment_beta_analysis_largersubset_list$fishless$randomized) %>% 
#   ggplot() +
#   geom_histogram(aes(x = qDTM.beta, fill = fish_info), 
#                  position="identity", alpha = 0.5, bins = 30) +
#   #scale_fill_manual(values = c('#005F73', '#AE2012')) +
#   scale_fill_manual(values = c('#0A9396', '#BB3E03')) +
#   #geom_vline(xintercept = within_treatment_beta_analysis_largersubset_list$fish$observed$qDTM.beta ) +
#   #geom_vline(xintercept = within_treatment_beta_analysis_largersubset_list$fishless$observed$qDTM.beta ) +
#   theme_bw() +
#   ylab('Count') + xlab('qDTM (beta)') +
#   theme(legend.title = element_blank())
# 
# ggsave(filename = here('figures', 'fish_fishless_beta_bootstrap.png'), 
#        width = 13*0.65, height = 9*0.65, bg = 'white')
# 
# 
# 
# 
# 
# rbind(within_treatment_beta_analysis_largersubset_list$fish$randomized,
#       within_treatment_beta_analysis_largersubset_list$fishless$randomized) %>% 
#   ggplot() +
#   geom_histogram(aes(x = M.beta.prime, fill = fish_info), 
#                  position="identity", alpha = 0.5, bins = 30) +
#   theme_bw()
# 
# 





# count_dat_processed_beta <- count_dat_processed %>% 
#   column_to_rownames(var = 'taxa') %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = 'lake') %>%
#   mutate(lake = gsub(" ", "_", str_trim(tolower(lake)))) %>%
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
# count_dat_processed_beta_fish <- count_dat_processed_beta %>% 
# filter(lake %in% processed_lake_info[processed_lake_info$fish == '1',]$lake) %>% 
#   column_to_rownames(var = 'lake')
# 
# count_dat_processed_beta_fishless <- count_dat_processed_beta %>% 
#   filter(lake %in% processed_lake_info[processed_lake_info$fish == '0',]$lake) %>% 
#   column_to_rownames(var = 'lake')



# count_dat_processed_beta_subset1 <- count_dat_processed_beta %>% 
#   filter(lake %in% c(pairwise_lake_combos[1, 1], pairwise_lake_combos[1, 2]) ) %>% 
#   column_to_rownames(var = 'lake')
# 
# count_dat_processed_beta_subset2 <- count_dat_processed_beta_subset1[,match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(count_dat_processed_beta_subset1) )]
# 
# count_dat_processed_beta$lake[!count_dat_processed_beta$lake %in% processed_lake_info$lake]
# 
# 
# beta_calc <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                       spmat = count_dat_processed_beta_subset2, 
#                       abund = FALSE, q = 1)




# rando_subsets_list <- replicate(5, count_dat_processed_beta %>%
#                                   column_to_rownames(var = 'lake') %>% 
#                                   slice_sample(n = 40, replace = FALSE),
#                                 simplify = FALSE)
# 
# 
# beta_rand_progress <- txtProgressBar(min = 0, max = length(rando_subsets_list), 
#                                      char = "*", style = 3)
# rand_beta_list <- list()
# 
# for (i in seq_along(rando_subsets_list)) {
#   
#   count_rand_subset <- rando_subsets_list[[i]][, match(rownames(as.matrix(zoo_dist_processed_func_space)), colnames(rando_subsets_list[[i]]) )]
#   
#   beta_calc_rand_subset <- FTD.beta(tdmat = zoo_dist_processed_func_space,
#                                     spmat = count_rand_subset,
#                                     abund = TRUE, q = 1)
#   
#   rand_beta_list[[paste('rand', i, sep = '_')]] <- as.data.frame((beta_calc_rand_subset[names(beta_calc_rand_subset) != 'disp.mat.weight']))
#   setTxtProgressBar(beta_rand_progress, i)
#   
# }
# 
# bind_rows(rand_beta_list, .id = 'rand_dataset')

# pairwise_beta_df %>% 
#   ggplot(aes(x = comparison_type, y = qDTM.beta)) +
#   geom_violin() +
#   geom_point(size = 0.5, alpha = 0.5, color = 'gray',
#              position = position_jitter(width = 0.1)) +
#   geom_boxplot(outlier.shape = NA, width = 0.1) +
#   theme_bw()
# 
# pairwise_beta_df %>% 
#   ggplot() +
#   geom_density(aes(x = qDTM.beta, color = comparison_type))

