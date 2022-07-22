######################################
######################################
### PROCESSING OF ZOOPLANKTON DATA ###
######################################
######################################


##########################
### SCRIPT PREPARATION ###
##########################

library(here)
library(tidyverse)
library(readxl)

### load data ###
#the rest of the data (length, counts) are loaded and processed together further down in the script
#egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length.xlsx'), sheet = 'Data'))
#egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_2_25_2022.xlsx'), sheet = 'Data'))
egg_len_dat <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_3_29_2022.xlsx'), sheet = 'Data formatted for R'))

rotifer_traits <- as.data.frame(read_excel(here('data', 'Rotifer_Measurements_and_Traits.xlsx'), na = "NA"))
trait_info_init <- as.data.frame(read_excel(here('data', 'Zooplankton_eggs_and_length_5_9_2022.xlsx'), sheet = 'Inter Traits'))

#COUNT DATA QUESTIONS
# #keratella, keratella 2 spine, keratella 1 spine in 2018
# #should d. middenorfiana/pulex and d. middenorfiana be treated separately or be lumped together?
# #difference between immature and juveline?
# 
# #Cyclapoid vs Cyclopoid unident. --> from the same lake


#######################
### PROCESSING DATA ###
#######################

### new processing step on 3/29/2022 ###
# rename to columns to convenient names ###
egg_len_dat_processed <- egg_len_dat %>% 
  rename("size_mm" = "Size (mm)",
         "egg_count" = "Number of eggs") %>% 
  rename_with(.fn = tolower, .cols = everything()) %>% 
  mutate(contains_fish = ifelse(fish == 1, 'yes', 'no'))


### 3/29/2022: THIS IS NOT LONGER NECESSARY WITH THE DATA REFORMATTING
### Length at reproduction and egg count data ###
#collected by Spencer Cruz
#Spencer is going to restructure the data when he is collecting it so that it is in "tidy" format, 
#so much of this processing may no longer be required for the egg count and length at reproduction data
# egg_len_dat$id <- seq_len(nrow(egg_len_dat)) #add unique ID for each row
# egg_len_dat$L_minutus_mm[egg_len_dat$L_minutus_mm == 's' & !is.na(egg_len_dat$L_minutus_mm)] <- NA
# egg_len_dat$L_minutus_mm <- as.numeric(egg_len_dat$L_minutus_mm)
# 
# egg_len_dat_processed <- egg_len_dat %>%
#   pivot_longer(cols = B_longirostris_mm:UnIdCyclopoid_eggs, names_to = 'species_trait', values_to = 'value', values_drop_na = TRUE) %>% 
#   mutate(species = gsub("_[A-Za-z]*$", '', species_trait),
#          trait = gsub("^[A-Za-z]*_*[A-Za-z]*_", '', species_trait),
#          Fish = if_else(fish == 1, 'fish', 'fishless')) %>% 
#   select(!c(species_trait, fish)) %>% 
#   pivot_wider(names_from = trait, values_from = value) %>% 
#   rename(reproduction_length = mm,
#          egg_count = eggs)


### LENGTH DATA ###
#collected by Lindsey Boyle
#length is of a collection of individuals of each taxa in the lake
#either all individuals were measured if only a few were in the sample or a random subset of individuals
#were measured if many individuals were collected (too many to be comprehensively measured)
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





trait_info <- trait_info_init %>% 
  select(Species, `Reproductive Mode`, `Body Shape`, `Feeding Type`) %>% 
  rename(species = Species,
         reproductive_mode = `Reproductive Mode`,
         body_shape = `Body Shape`,
         feeding_type = `Feeding Type`)

length_dat_processed <- zoop_length_total_list0 %>%
  bind_rows() %>% 
  group_by(taxa) %>% 
  summarise(mean_length = mean(length_mm, na.rm = TRUE), .groups = 'drop')

trait_info <- trait_info %>% 
  na.omit() %>% 
  mutate(taxa = case_when(species == 'B_longirostris' ~ 'b. longirostris',
                          species == 'D_mendotae' ~ 'd. mendatoe',
                          species == 'D_pulex' ~ 'd. middenorfiana/pulex',
                          species == 'H_shoshone' ~ "h. shoshone",
                          species == 'L_minutus' ~ 'l. minutus',
                          species == 'Chydorus' ~ 'chydorus',
                          species == 'H. gibberum' ~ 'h. gibberum',
                          species == 'P. pediculous' ~ 'p. pediculus')) %>%
  select(!species) %>% 
  left_join(., length_dat_processed, by = 'taxa') %>% 
  rename(species = taxa)
  #column_to_rownames(var = "taxa")


rotifer_info <- rotifer_traits %>% 
  select(Species, `Reproductive Mode`, `Body Shape`, `Trophic Group`, `Feeding Type`, `Size(mm)`) %>% 
  rename(species = Species,
         reproductive_mode = `Reproductive Mode`,
         body_shape = `Body Shape`,
         feeding_type = `Feeding Type`,
         length = `Size(mm)`) %>% 
  group_by(species) %>% 
  summarize(reproductive_mode = first(reproductive_mode),
            body_shape = first(body_shape),
            feeding_type = first(feeding_type),
            mean_length = mean(length)) %>% 
  mutate(species = tolower(species)) %>% 
  mutate(species = case_when(species == "colonial_conochilidae" ~ "colonial conochilidae",
                             species == "fillinia" ~ "filinia",
                             species == "kellicotia" ~ "kellicotta",
                             species == 'keratella_2spine' ~ "keratella 2 spine",
                             species == "keratella_1spine" ~ "keratella 1 spine",
                             species == "synchaeta" ~ "syncheata",
                             species == "asplancha" ~ "asplancha",
                             species == "notholca" ~ "notholca",
                             species == "polyarthra" ~ "polyarthra",
                             species == "euchlanis" ~ "euchlanis"))


all_taxa_trait_info <- rbind(trait_info, rotifer_info)

all_taxa_trait_info$species[!(all_taxa_trait_info$species %in% unique(c(count_dat_processed_list$length_dat_2018$taxa, count_dat_processed_list$length_dat_2019$taxa) ))]

# count_dat_processed_list_name_update <- lapply(count_dat_processed_list, function(x) {
#   x %>% 
#     mutate(Lake = recode(Lake, "north_blue" = "n of blue"))
# })


### COUNT DATA ###
#collected by Lindsey Boyle
#The number of individuals for each taxa. Some of the taxa have more specific count information such as
#counts by sex, counts of individuals with eggs, etc...
count_dat_processed_list <- list()

for (YEAR in c(2018, 2019)) {
  
  count_dat_processed_list_year <- lapply(setNames(nm = excel_sheets(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')))), function(NAME) {
    message(NAME)
    count_dat <- as.data.frame(read_excel(here('data', paste0('WR', YEAR, 'ZoopsID.xlsx')), sheet = NAME, skip = 1, n_max = 1, na = 'NA')) #read in the dataframe (skip first 5 rows, which has count info)
    count_dat_1 <- count_dat[, !grepl(pattern = "^\\.{3}[0-9]*", colnames(count_dat))] #remove any extra columns (named "...NUMBER")
    
    count_dat_2 <- count_dat_1 %>% 
      pivot_longer(cols = everything(), names_to = 'taxa_counttype', values_to = 'count') %>% #pivot data so that count and taxa are in columns
      mutate(taxa_counttype = str_replace(taxa_counttype, pattern = 'juvenile|juveline', replace = 'immature'), #switch juvenile and a mispelling of that to immature
             count_type = str_trim(tolower(str_extract(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|w/[ ]*eggs$|w/[ ]*epiphia$|immature$"))), #add count type (e.g. m, f, total) to new column
             taxa = str_trim(tolower(str_remove(taxa_counttype, pattern = "[Tt]otal$|[Ff]$|[Mm]$|[Ff]*[ ]*w/[ ]*eggs$|[Ff]*[ ]*w/[ ]*epiphia$|immature$|\\?")))) %>% #add taxa to new column
      mutate(taxa = str_remove(taxa, '\\.$|\\?'), #remove trailing periods and question marks from taxa name
             count_type = str_remove(count_type, pattern = '\\s*w/\\s*'), #remove "w/" from count type (allowing for matching with differences in spaces)
             count_type = recode(count_type, juvenile = 'immature', juveline = 'immature'),
             count_type = if_else(is.na(count_type), 'total', count_type)) %>% #add total as count type for rows that don't yet have count type
      mutate(taxa = if_else(str_detect(string = taxa, pattern = 'd\\. men*dotae'), 'd. mendotae', taxa),
             taxa = if_else(taxa == "d. mendotae", "d. mendatoe", taxa),
             taxa = if_else(taxa == "fillinia", "filinia", taxa),
             taxa = if_else(taxa %in% c("euclanis sp", "euchlanis sp"), "euchlanis", taxa),
             taxa = if_else(taxa == "h. gibberu", "h. gibberum", taxa),
             taxa = if_else(taxa == "d. middenorfiana", "d. middenorfiana/pulex", taxa),
             taxa = ifelse(taxa %in% c("polyphemus pediculus"), "p. pediculus", taxa)) %>%
      select(!taxa_counttype) %>% 
      pivot_wider(names_from = count_type, values_from = count, values_fill = 0) %>%
      rowwise() %>% 
      mutate(total = ifelse(total == 0, ifelse(m != 0 | f != 0, sum(m, f, na.rm = TRUE), sum(eggs, na.rm = TRUE)), total)) %>% 
      ungroup() %>% 
      as.data.frame()
    
    #add immature column if doesnt exist (only one year identified immature/juvenile)
    if (!('immature' %in% colnames(count_dat_2)))
      count_dat_2 <- count_dat_2 %>% add_column(immature = 0, .after = "epiphia")
    
    #replace NAs in the count columns with 0 and remove rows where all counts are 0
    count_dat_3 <- count_dat_2 %>% 
      replace_na(list(total = 0, 
                      m = 0, 
                      f = 0, 
                      eggs = 0, 
                      immature = 0)
      ) %>%
      filter(if_any(c(total, m, f, eggs, immature), ~ . != 0))
    
    return(count_dat_3)
  })
  
  #unique(sapply(count_dat_2018_processed_list, ncol))
  
  count_dat_processed_list[[paste0('length_dat_', YEAR)]] <- bind_rows(count_dat_processed_list_year, .id = 'Lake')
  count_dat_processed_list[[paste0('length_dat_', YEAR)]]$year <- YEAR
}
# prac_df <- data.frame(total = c(1, 0, 1, 0, 0, 2),
#            a     = c(10, 2, 1, 0, 0, 2),
#            b     = c(10, 3, 1, 0, 0, 2),
#            c     = c(10, 2, 1, 100, NA, 2))
# 
# prac_df %>%
#   rowwise() %>% 
#   mutate(total = ifelse(total == 0, ifelse(a != 0 | b != 0, sum(a, b, na.rm = TRUE), sum(c, na.rm = TRUE)), total)) %>% 
#   ungroup()

saveRDS(all_taxa_trait_info, here('data', 'processed_data', 'all_taxa_trait_info_processed.rds') )
saveRDS(count_dat_processed_list, here('data', 'processed_data', 'count_dat_processed_list.rds') )
saveRDS(zoop_length_total_list0, here('data', 'processed_data', 'zoop_length_total_list0.rds') )



