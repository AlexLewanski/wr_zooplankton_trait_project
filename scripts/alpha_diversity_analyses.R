####################################
####################################
### ANALYSIS OF ZOOPLANKTON DATA ###
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


### QUESTIONS ###
#difference in variables across years
#include depth?
#NDVI?
#species-specific effects?



##########################
### SCRIPT PREPARATION ###
##########################

### load libraries ###
library(here)
library(tidyverse)
library(gawdis)
library(readxl)
library(picante) #matrix randomization
#library(PerformanceAnalytics)
#library(ggcorrplot)
library(MuMIn)
library(GGally) #correlation plot
library(lubridate) #handling times and dates
library(glmmTMB) #mixed model package
library(AICcmodavg) #model comparison
library(performance) #used for vif calc (alternative is car package but vif func. doesn't work with glmmTMB)
library(gstat) #semivariogram
library(ape) #variogram
library(DHARMa) #diangostics for mixed models
#library(AER) #dispersion test for poisson glm


### custom functions
source(here("scripts", "custom_functions_zooplankton_project.R"))


### loading data (and initial processing) ###
lindsey_master_data <- as.data.frame(read_excel(here('data', 'WR_data_ALL_Boyle.xlsx')))
fish_info <- lindsey_master_data #%>% 
  #select(lake, fish)

processed_lake_info <- fish_info[!duplicated(fish_info$lake),] %>%
  #mutate(long = abs(long)*-1)
   mutate(long = abs(long)*-1,
          lake = dplyr::recode(lake,
                        'lake37' = "lake_37",
                        'blackjoe' = "black_joe",
                        'upper_blackjoe' = "upper_black_joe",
                        "below_tb" = "lake_below_tb"))

alpha_ses_fd <- read.csv(file = here('results', 'alpha_ses_fd.csv'))
alpha_scheiner_fd <- read.csv(file = here('results', 'alpha_scheiner_metrics.csv'))


alpha_scheiner_fd_lake <- left_join(alpha_scheiner_fd, processed_lake_info, by = 'lake') %>%
  select(!q) %>% 
  #filter(lake != 'lost_lake') %>% 
  #select(!community) %>% 
  pivot_longer(cols = nsp:qDTM,
               names_to = "fd_metric", 
               values_to = "value") %>%
    select(lake, lat, long, elev, area, #trip, region, side, 
           date, time, trip, conditions, temp, ph, 
           cond, do, secci, maxd, tow, vol_sampled,
           fish, bk, ct, gd, ctgd, rb, splk,
           P_mg, N_mg, C_mg, cla, P_cla, DOC,
           Lithium, Sodium, Ammonium, Potassium, 
           Magnesium, Calcium, Fluoride, Chloride, 
           Nitrite, Nitrate, Phosphate, 
           fd_metric, value) %>% 
    mutate(trip = as.character(trip),
           fish = as.character(fish),
           julian_date = yday(date),
           time_convert = hour(time) + (minute(time)/60) + (second(time)/3600))


fd_ses_df_lake <- left_join(alpha_ses_fd, processed_lake_info, by = 'lake')  %>%
  #filter(lake != 'lost_lake') %>% 
  #select(!community) %>% 
  pivot_longer(cols = ends_with('ses'),
               names_to = "fd_metric", 
               values_to = "value") %>%
  select(lake, lat, long, elev, area, #trip, region, side, 
         date, time, trip, conditions, temp, ph, 
         cond, do, secci, maxd, tow, vol_sampled,
         fish, bk, ct, gd, ctgd, rb, splk,
         P_mg, N_mg, C_mg, cla, P_cla, DOC,
         Lithium, Sodium, Ammonium, Potassium, 
         Magnesium, Calcium, Fluoride, Chloride, 
         Nitrite, Nitrate, Phosphate, 
         fd_metric, value) %>% 
  mutate(trip = as.character(trip),
         fish = as.character(fish),
         julian_date = yday(date),
         time_convert = hour(time) + (minute(time)/60) + (second(time)/3600))

# fd_ses_df_lake_wider <- fd_ses_df_lake %>% 
#   pivot_wider(names_from = fd_metric, values_from = value) %>% 
#   as.data.frame()

fd_combined <- rbind(fd_ses_df_lake, alpha_scheiner_fd_lake)

fd_ses_df_lake_wider <- fd_combined %>% 
  pivot_wider(names_from = fd_metric, values_from = value) %>% 
  mutate(temp_std = standardize_val(temp),
         area_std = standardize_val(area),
         time_convert_std = standardize_val(time_convert),
         julian_date_std = standardize_val(julian_date)) %>%  
  as.data.frame()



##############################################
### Step 1: Are there outliers in Y and X? ###
##############################################

#boxplot
fd_ses_df_lake_wider %>% 
  select(lake, fric_ses, fdis_ses, fdiv_ses, feve_ses, julian_date, time_convert, area, elev, M.prime, qEt, nsp, qDTM) %>% 
  pivot_longer(cols = c(fric_ses, fdis_ses, fdiv_ses, feve_ses, julian_date, time_convert, area, elev, M.prime, qEt, nsp, qDTM), 
               names_to = 'variable', values_to = 'value') %>%
  ggplot() +
  geom_boxplot(aes(y = value)) +
  facet_wrap(~ variable, scales = 'free') +
  theme_bw()


#cleveland dotplot
fd_ses_df_lake_wider %>% 
  select(lake, fric_ses, fdis_ses, fdiv_ses, feve_ses, julian_date, time_convert, area, elev, M.prime, qEt, nsp, qDTM) %>% 
  mutate(row_num = 1:n()) %>% 
  pivot_longer(cols = c(fric_ses, fdis_ses, fdiv_ses, feve_ses, julian_date, time_convert, area, elev, M.prime, qEt, nsp, qDTM), 
               names_to = 'variable', values_to = 'value') %>%
  ggplot() +
  geom_point(aes(x = value, y = row_num)) +
  facet_wrap(~ variable, scales = 'free') +
  theme_bw()



###################################################
### Step 2: Do we have homogeneity of variance? ###
###################################################

var_text <- fd_combined %>% 
  filter(fd_metric %in% c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM')) %>% 
  mutate(fd_metric = factor(fd_metric, levels = c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM'))) %>% 
  group_by(fish, fd_metric) %>% 
  summarize(variance = var(value), .groups = 'drop') %>% 
  group_by(fd_metric) %>% 
  summarize(var_ratio = round(max(variance)/min(variance), 2) )

fd_combined %>% 
  filter(fd_metric %in% c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM')) %>% 
  mutate(fd_metric = factor(fd_metric, levels = c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM'))) %>% 
  ggplot() +
  geom_boxplot(aes(x = fish, y = value)) +
  facet_wrap(~fd_metric, scales = 'free') +
  theme_bw() +
  geom_text(data = var_text %>% mutate(label = paste0('var. ratio = ', var_ratio)), 
            mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust   = -0.1,
            vjust   = -1)



##################################################
### Step 3: Are the data normally distributed? ###
##################################################

fd_ses_df_lake_wider %>% 
  select(lake, fric_ses, fdis_ses, fdiv_ses, feve_ses, M.prime, qEt, nsp, qDTM) %>% 
  mutate(row_num = 1:n() ) %>% 
  pivot_longer(cols = c(fric_ses, fdis_ses, fdiv_ses, feve_ses, M.prime, qEt, nsp, qDTM), 
               names_to = 'variable', values_to = 'value') %>%
  ggplot() +
  geom_histogram(aes(x = value), bins = 15) +
  facet_wrap(~ variable, scales = 'free') +
  theme_bw()



###########################################################
### Step 5: Is there collinearity among the covariates? ###
###########################################################

#https://stackoverflow.com/questions/61686171/how-to-add-the-spearman-correlation-p-value-along-with-correlation-coefficient-t

#ggpairs(fd_ses_df_lake_wider[c('temp', 'area', 'julian_date', 'time_convert', 'elev')], 
#        title="correlogram of predictor variables",
#        upper = list(continuous = wrap("cor", method = "spearman"))) +
#  theme_bw()

fd_ses_df_lake_wider %>% 
  select(temp, area, julian_date, time_convert, elev, fish) %>%
  mutate(fish = as.numeric(fish))

fd_ses_df_lake_wider %>% 
  select(temp, area, julian_date, time_convert, elev, fish) %>%
  mutate(fish = as.numeric(fish)) %>% 
  ggpairs(.,
       title="correlogram of predictor variables",
       upper = list(continuous = wrap("cor", method = "pearson"))) +
 theme_bw()

#Point-biserial correlation
#cor.test(as.numeric(fd_ses_df_lake_wider$fish), fd_ses_df_lake_wider$elev)

fd_ses_df_lake_wider %>% 
  select(temp, area, julian_date, time_convert, elev, fish) %>%
  mutate(fish = as.numeric(fish)) %>% 
  cor.test(.x$fish, .x$temp)


yos <- sample(1:30, 100, replace = TRUE)
disease <- sample(0:1, 100, replace = TRUE)
cor.test(yos, disease)



#####################################################################
### Step 6: What are the relationships between Y and X variables? ###
#####################################################################

scatterplot_list <- lapply(setNames(nm = c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM')), function(fd, dataset) {
  fd_subset <- gsub("_ses", "", fd)
  fd_title <- paste0(toupper(substr(fd_subset, 1, 2)), substr(fd_subset, 3, nchar(fd_subset)))
  
  return(
    dataset %>% 
      select(lake, fd_metric, value, temp, area, julian_date, time_convert, fish) %>% 
      mutate(fish = as.numeric(fish)) %>% 
      pivot_longer(cols = temp:fish, names_to = 'predictor', values_to = 'predictor_val') %>% 
      filter(fd_metric == fd) %>% 
      ggplot() +
      geom_point(aes(predictor_val, value)) +
      facet_wrap(~predictor, scales = 'free') +
      ggtitle(fd_title) + 
      xlab('Predictor') +
      ylab(fd_title) +
      theme_bw()
  )
  
}, dataset = fd_combined)

scatterplot_list$fric_ses
scatterplot_list$fdis_ses
scatterplot_list$fdiv_ses
scatterplot_list$feve_ses
scatterplot_list$M.prime
scatterplot_list$qEt
scatterplot_list$nsp
scatterplot_list$qDTM



###################################
### ADDITIONAL DATA EXPLORATION ###
###################################

### how does date correspond to trip identity ###
fd_ses_df_lake_wider %>% 
  ggplot() +
  geom_histogram(aes(x = julian_date, fill = trip), 
                 position = "identity", alpha = 0.7) +
  theme_bw()



############################################################################
### MODEL FITTING 1: (generalized) linear models with date as predictor  ###
############################################################################

### STEP 1. FIT MODELS ####
predictor_vec <- setNames(c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM'), 
                          nm = c('FRic', 'FDis', 'FDiv', 'FEve', 'M.prime', 'qEt', 'nsp', 'qDTM'))

alpha_fd_mods <- lapply(predictor_vec, function(outcome, data) {
  
  family <- switch(outcome, nsp = {'poisson'}, 'gaussian')
  
  fd_global_mod <- glm(as.formula(paste0(outcome, ' ~ temp + area + julian_date + time_convert + fish')),
                       data = fd_ses_df_lake_wider, family = family, na.action = "na.fail")
  
  #fd_global_mod <- lm(as.formula(paste0(outcome, ' ~ temp + area + julian_date + time_convert + fish')),
  #                    data = fd_ses_df_lake_wider, na.action = "na.fail")
                  
  mod_dredge <- fric_dredge <- dredge(global.model = fd_global_mod, 
                                      beta = c("none"), 
                                      evaluate = TRUE,
                                      rank = "AICc")
  
  mod_extract <- get.models(fric_dredge, subset = TRUE) #subset = delta <= 3
  
  return(list(global_mod = fd_global_mod,
              dredge_output = mod_dredge,
              top_mods = mod_extract) 
         )
  
}, data = fd_ses_df_lake_wider)

mod_rename_list_glm <- list()
for (FD in names(alpha_fd_mods)) {
  for (x in names(alpha_fd_mods[[FD]][["top_mods"]])) {
    mod_var <- paste(rownames(summary(alpha_fd_mods[[FD]][["top_mods"]][[x]])$coefficients), collapse = " + ")
    mod_rename <- gsub(pattern = '\\(Intercept\\)', 'intercept', mod_var)
    mod_rename_list_glm[[FD]][[mod_rename]] <- alpha_fd_mods[[FD]][["top_mods"]][[x]]
  }
}

selection_table_list <- lapply(mod_rename_list_glm, function(x) aictab(cand.set = x) )



######################################
### POST-MODEL EVALUATION --- GLMS ###
#####################################

### 1. CHECKING FOR PROBLEMATIC MULTICOLLINEARITY ###
#https://stats.stackexchange.com/questions/10079/rules-of-thumb-for-minimum-sample-size-for-multiple-regression
glm_vif_list <- lapply(setNames(nm = c('FRic', 'FDis', 'FDiv', 'FEve', 'M.prime', 'qEt', 'nsp', 'qDTM')), function(x, mod_list) {
  #vif(mod_list[[x]]$global_mod)
  check_collinearity(mod_list$FRic$global_mod)
}, mod_list = alpha_fd_mods)


### 2. MODEL DIAGNOSTICS ###
glm_diagnostics(mod_list = mod_rename_list_glm,
                file_name = here('results', 'mod_diagnostics', 'glm_diagnostic_plots'),
                panel_width = 8,
                plot_height = 8,
                show_progress = TRUE)

#https://stats.stackexchange.com/questions/70558/diagnostic-plots-for-count-regression
#https://stats.stackexchange.com/questions/66586/is-there-a-test-to-determine-whether-glm-overdispersion-is-significant/66593#66593
glm_dispersion_test_df <- lapply(setNames(nm = names(mod_rename_list_glm$nsp)), function(x, mod_list) {
  dispersion_test <- dispersiontest(mod_list[[x]])
  
  return(data.frame(statistic = glm_dispersion_list$`intercept + area`$statistic,
                    p = glm_dispersion_list$`intercept + area`$p.value,
                    estimate = glm_dispersion_list$`intercept + area`$estimate,
                    null_value = glm_dispersion_list$`intercept + area`$null.value,
                    alternative = glm_dispersion_list$`intercept + area`$alternative,
                    method = glm_dispersion_list$`intercept + area`$method)
         )
}, mod_list = mod_rename_list_glm$nsp) %>% 
  bind_rows(.id = 'model')


### 3. CHECKING FOR EVIDENCE OF SPATIAL AUTOCOR IN EACH MODEL ###
glm_sp_autocor_list <- list()
for (i in c('moran', 'variogram')) {
  glm_sp_autocor_list[[i]] <- lapply(mod_rename_list_mixed, function(FD, data, i) {
    lapply(FD, function(model_object, data, i) {
      eval_spatial_autocor(data = data, 
                           mod = model_object, 
                           long_col = 'long', 
                           lat_col = 'lat', 
                           method = i)
    }, data = data, i = i)
  }, data = fd_ses_df_lake_wider, i = i)
}

#multipanel of semivariograms for each outcome variable's set of models
glm_variogram_multipan_list <- lapply(setNames(nm = names(glm_sp_autocor_list$variogram)), function(x, var_results) {
  lapply(var_results[[x]], function(y) y[['variogram']]) %>% 
    bind_rows(.id = 'model') %>% 
    ggplot() +
    geom_point(aes(x = dist, y = gamma)) + 
    facet_wrap(~model) +
    theme_bw()
}, var_results = glm_sp_autocor_list$variogram)

#combined dataframe of all morgan test results
glm_moran_combined_df <- lapply(setNames(nm = names(glm_sp_autocor_list$moran)), function(FD, moran_list) {
  lapply(moran_list[[FD]], function(x) as.data.frame(x$moran)) %>% 
    bind_rows(.id = 'predictors')
}, moran_list = glm_sp_autocor_list$moran) %>% 
  bind_rows(.id = 'outcome') %>% 
  mutate(model = paste(outcome, '~', predictors, sep = " "))



############################################################
### MODEL FITTING 2: (g)lmms with trip as random effect  ###
############################################################

predictor_vec <- setNames(c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM'), 
                          nm = c('FRic', 'FDis', 'FDiv', 'FEve', 'M.prime', 'qEt', 'nsp', 'qDTM'))


### STEP 1. TEST IF POISSON DISTRIBUTION WORKS FOR THE SPECIES RICHNESS MODELS ####

nsp_fd_mods_mixed_poisson <- lapply(c(nsp = 'nsp'), function(outcome, data) {

  fd_global_mod <- glmmTMB(as.formula(paste0(outcome, ' ~ temp_std + area_std + time_convert_std + fish + (1|trip)')),
                           data = data, family = 'poisson', na.action = "na.fail")
  
  mod_dredge <- dredge(global.model = fd_global_mod, 
                       beta = c("none"), 
                       evaluate = TRUE,
                       rank = "AICc")
  
  mod_extract <- get.models(mod_dredge, subset = TRUE) #subset = delta <= 3
  
  return(list(global_mod = fd_global_mod,
              dredge_output = mod_dredge,
              top_mods = mod_extract) 
  )
  
}, data = fd_ses_df_lake_wider)

#is there evidence of over/underdispersion based on the poisson glmm of species richness?
#answer: YES, there is evidence of underdispersion in all models
#let's switch to conway-maxwell-poisson distribution
test_disp_list_poisson <- lapply(nsp_fd_mods_mixed_poisson$nsp$top_mods, function(x) {
  testDispersion(simulateResiduals(x))
})


### STEP 2. FIT MODELS ####
alpha_fd_mods_mixed <- lapply(predictor_vec, function(outcome, data) {
  message(outcome)
  
  family <- switch(outcome, nsp = {'compois'}, 'gaussian')
  
  #fd_global_mod <- glmer(as.formula(paste0(outcome, ' ~ temp + area + time_convert + fish + (1|trip)')),
  #                     data = fd_ses_df_lake_wider, family = family, na.action = "na.fail")
  
  fd_global_mod <- glmmTMB(as.formula(paste0(outcome, ' ~ temp_std + area_std + time_convert_std + fish + (1|trip)')),
                         data = data, family = family, na.action = "na.fail")
  
  mod_dredge <- dredge(global.model = fd_global_mod, 
                                      beta = c("none"), 
                                      evaluate = TRUE,
                                      rank = "AICc")
  
  mod_extract <- get.models(mod_dredge, subset = TRUE) #subset = delta <= 3
  
  return(list(global_mod = fd_global_mod,
              dredge_output = mod_dredge,
              top_mods = mod_extract) 
  )
  
}, data = fd_ses_df_lake_wider)


mod_rename_list_mixed <- list()
for (FD in names(alpha_fd_mods_mixed)) {
  for (x in names(alpha_fd_mods_mixed[[FD]][["top_mods"]])) {
    mod_var <- paste(rownames(summary(alpha_fd_mods_mixed[[FD]][["top_mods"]][[x]])$coefficients$cond), collapse = " + ")
    mod_rename <- gsub(pattern = '\\(Intercept\\)', 'intercept', mod_var)
    mod_rename_list_mixed[[FD]][[mod_rename]] <- alpha_fd_mods_mixed[[FD]][["top_mods"]][[x]]
  }
}

selection_table_list_mixed <- lapply(mod_rename_list_mixed, function(x) aictab(cand.set = x) )


#does the conway-maxwell-poisson distribution deal with the underdispersion present in the poisson models?
test_disp_list_compois <- lapply(alpha_fd_mods_mixed$nsp$top_mods, function(x) {
  testDispersion(simulateResiduals(x))
})



##########################################
### GLMMS --- ONLY FISH AS A PREDICTOR ###
##########################################

predictor_vec <- setNames(c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses', 'M.prime', 'qEt', 'nsp', 'qDTM'), 
                          nm = c('FRic', 'FDis', 'FDiv', 'FEve', 'M.prime', 'qEt', 'nsp', 'qDTM'))

### STEP 1. FIT MODELS ####

alpha_fd_mods_mixed_fish <- lapply(predictor_vec, function(outcome, data) {
  message(outcome)
  
  #family <- switch(outcome, nsp = {'poisson'}, 'gaussian')
  family <- switch(outcome, nsp = {'compois'}, 'gaussian')
  
  #fd_global_mod <- glmer(as.formula(paste0(outcome, ' ~ temp + area + time_convert + fish + (1|trip)')),
  #                     data = fd_ses_df_lake_wider, family = family, na.action = "na.fail")
  
  #fd_fish_mod <- glmmTMB(as.formula(paste0(outcome, ' ~ fish + (1|trip)')),
  #                         data = data, family = family, na.action = "na.fail")
  fd_fish_mod <- glmmTMB(as.formula(paste0(outcome, ' ~ fish + (1|trip)')),
                         data = data, family = family, na.action = "na.fail")
  
  return(fd_fish_mod)
  
}, data = fd_ses_df_lake_wider)


### VISUALIZING RESULTS ###
#NS: P > 0.10,
#.: P < 0.10,
#*: P < 0.05,
#**: P < 0.01,
#***: P < 0.001,

fish_mod_summary <- lapply(setNames(nm = names(alpha_fd_mods_mixed_fish)), function(x, mod_list) {

  return(
    data.frame(estimate = summary(mod_list[[x]])$coefficients$cond['fish1', 'Estimate'],
               p = summary(mod_list[[x]])$coefficients$cond['fish1', 'Pr(>|z|)']) %>%
      mutate(estimate_sign = if_else(estimate > 0, '+', '-'),
             p_symbol = case_when(p > 0.10 ~ 'ns',
                                  p <= 0.10 & p > 0.05 ~ '\U2022', #interpunct instead of period (.)
                                  p <= 0.05 & p > 0.01 ~ '*',
                                  p <= 0.01 & p > 0.001 ~ '**',
                                  p <= 0.001 ~ '**'),
             #figure_text = paste0(p_symbol, ' (', estimate_sign, ')'),
             color = if_else(estimate > 0, 'green', 'red'),
             figure_text = paste0(estimate_sign, '\n', p_symbol)
             )
  )

}, mod_list = alpha_fd_mods_mixed_fish) %>% 
  bind_rows(.id = 'metric')


fd_longer_format <- fd_ses_df_lake_wider %>% 
  select(lake, fish, fric_ses, fdis_ses, fdiv_ses, feve_ses, nsp, M.prime, qEt, qDTM) %>% 
  pivot_longer(cols = fric_ses:qDTM, names_to = 'metric', values_to = 'value')

multidim_fd_alpha_plots <- fd_longer_format %>%
  filter(metric %in% c('fric_ses', 'fdis_ses', 'fdiv_ses', 'feve_ses')) %>% 
  mutate(fish = if_else(fish == '1', 'Fish', 'Fishless'),
         metric = factor(dplyr::recode(metric,
                                'fric_ses' = 'Richness',
                                'fdis_ses' = 'Dispersion',
                                'fdiv_ses' = 'Divergence',
                                'feve_ses' = 'Evenness'),
                         levels = c('Richness', 'Dispersion', 'Divergence', 'Evenness')),
         fish = factor(fish, levels = c('Fishless', 'Fish'))) %>%
  ggplot() +
  #geom_boxplot(aes(x = metric, y = value, fill = fish),
  #             alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
  geom_boxplot(aes(x = metric, y = value, fill = fish, color = fish),
               alpha = 0.2, lwd = 0.8, outlier.shape = NA, #color = '#a6a6a6', 
               position = position_dodge(width = 0.875)) +
  #geom_violin(aes(x = metric, y = value, fill = fish),
  #            alpha = 0.3, color = '#a6a6a6', lwd = 0.8) +
  geom_point(aes(x = metric, y = value, fill = fish, group = fish, color = fish),
             size = 2,
             position = position_jitterdodge(dodge.width = 0.875, jitter.width = 0.2)) +
  geom_text(data = fish_mod_summary %>%
              filter(metric %in% c('FRic', 'FDis', 'FDiv', 'FEve')) %>% 
              mutate(metric = factor(dplyr::recode(metric,
                                                   'FRic' = 'Richness',
                                                   'FDis' = 'Dispersion',
                                                   'FDiv' = 'Divergence',
                                                   'FEve' = 'Evenness'),
                                     levels = c('Richness', 'Dispersion', 'Divergence', 'Evenness'))), 
            aes(x = metric, y = 2.6, label = figure_text), 
            col = 'black', size = 8, lineheight = 0.75)  +
  xlab('Functional diversity') +
  ylab('Standardized effect size') +
  scale_fill_manual(values = c('#DC3220', '#005AB5')) +
  scale_color_manual(values = c('#DC3220', '#005AB5')) +
  theme_bw() +
  theme(axis.title = element_text(size = 19),
        axis.text.x = element_text(size = 15),
        legend.title=element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  ylim(-2, 2.75)

ggsave(multidim_fd_alpha_plots,
       filename = here('figures', 'multidim_fd_alpha.png'), 
       width = 13*0.8, height = 7*0.8)
  

scheiner_alpha_plot <- fd_longer_format %>%
  filter(metric %in% c('nsp', 'M.prime', 'qEt', 'qDTM')) %>% 
  mutate(fish = factor(if_else(fish == '1', 'Fish', 'Fishless'), levels = c('Fishless', 'Fish')),
         metric = factor(dplyr::recode(metric,
                                       'nsp' = 'Species richness',
                                       'M.prime' = 'M.prime',
                                       'qEt' = 'qEt',
                                       'qDTM' = 'qDTM'),
                         levels = c('qDTM', 'M.prime', 'qEt', 'Species richness'))) %>%
  ggplot() +
  #geom_boxplot(aes(x = metric, y = value, fill = fish),
  #             alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8) +
  geom_boxplot(aes(x = metric, y = value, fill = fish),
               alpha = 0.3, outlier.shape = NA, color = '#a6a6a6', lwd = 0.8, position = position_dodge(width = 0.875)) +
  #geom_violin(aes(x = metric, y = value, fill = fish),
  #            alpha = 0.3, color = '#a6a6a6', lwd = 0.8) +
  geom_point(aes(x = metric, y = value, fill = fish, group = fish, color = fish),
             size = 2,
             position = position_jitterdodge(dodge.width = 0.875, jitter.width = 0.2)) +
  geom_text(data = fish_mod_summary %>%
              filter(metric %in% c('nsp', 'M.prime', 'qEt', 'qDTM')) %>% 
              mutate(metric = factor(dplyr::recode(metric,
                                                   'nsp' = 'Species richness',
                                                   'M.prime' = 'M.prime',
                                                   'qEt' = 'qEt',
                                                   'qDTM' = 'qDTM'),
                                     levels = c('qDTM', 'M.prime', 'qEt', 'Species richness'))), 
            aes(x=metric, y = Inf, label = figure_text), col = 'black', 
            size = 8, vjust = 1, lineheight = 0.75) +
  xlab('Functional diversity') +
  ylab('Standardized effect size') +
  scale_fill_manual(values = c('#DC3220', '#005AB5')) +
  scale_color_manual(values = c('#DC3220', '#005AB5')) +
  theme_bw() +
  theme(#axis.title = element_text(size = 19),
        axis.title = element_blank(),
        #axis.text.x = element_text(size = 15),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~metric, scales = 'free')

ggsave(scheiner_alpha_plot,
       filename = here('figures', 'scheiner_fd_alpha.png'), 
       width = 13*0.8, height = 7*0.8)


  

mod_rename_list_mixed <- list()
for (FD in names(alpha_fd_mods_mixed)) {
  for (x in names(alpha_fd_mods_mixed[[FD]][["top_mods"]])) {
    mod_var <- paste(rownames(summary(alpha_fd_mods_mixed[[FD]][["top_mods"]][[x]])$coefficients$cond), collapse = " + ")
    mod_rename <- gsub(pattern = '\\(Intercept\\)', 'intercept', mod_var)
    mod_rename_list_mixed[[FD]][[mod_rename]] <- alpha_fd_mods_mixed[[FD]][["top_mods"]][[x]]
  }
}

selection_table_list_mixed <- lapply(mod_rename_list_mixed, function(x) aictab(cand.set = x) )





#######################################
### POST-MODEL EVALUATION --- GLMMS ###
#######################################

### 1. CHECKING FOR PROBLEMATIC MULTICOLLINEARITY
#https://stats.stackexchange.com/questions/10079/rules-of-thumb-for-minimum-sample-size-for-multiple-regression
mixedmod_vif_list <- lapply(setNames(nm = c('FRic', 'FDis', 'FDiv', 'FEve', 'M.prime', 'qEt', 'nsp', 'qDTM')), function(x, mod_list) {
  check_collinearity(mod_list[[x]]$global_mod)
}, mod_list = alpha_fd_mods_mixed)


### 2. MODEL DIAGNOSTICS ###
mixed_mod_diagnostics(mod_list = mod_rename_list_mixed,
                      file_name = here('results', 'mod_diagnostics', 'mixed_mod_diagnostic_plots'),
                      diagnostics = c('QQ', 'resid_vs_predict', 'dispersion'),
                      panel_width = 10,
                      plot_height = 6,
                      show_progress = TRUE)


### 3. CHECKING FOR EVIDENCE OF SPATIAL AUTOCOR IN EACH MODEL ###
mixed_mod_sp_autocor_list <- list()
for (i in c('moran', 'variogram')) {
  
  mixed_mod_sp_autocor_list[[i]] <- lapply(mod_rename_list_mixed, function(FD, data, i) {
    lapply(FD, function(model_object, data, i) {
      switch(i, 
             variogram = {variogram_calc(data = data, mod = model_object, long_col = 'long', lat_col = 'lat')},
             moran = {moran_wrapper_updated(data = data, mod = model_object, long_col = 'long', lat_col = 'lat')})
    }, data = data, i = i)
  }, data = fd_ses_df_lake_wider, i = i)
}

#multipanel of semivariograms for each outcome variable's set of models
mixed_mod_variogram_multipan_list <- lapply(setNames(nm = names(mixed_mod_sp_autocor_list$variogram)), function(x, var_results) {
  var_results[[x]] %>% 
    bind_rows(.id = 'model') %>% 
    ggplot() +
    geom_point(aes(x = dist, y = gamma)) + 
    facet_wrap(~model) +
    theme_bw()
}, var_results = mixed_mod_sp_autocor_list$variogram)

#combined dataframe of all morgan test results
mixed_mod_moran_combined_df <- lapply(setNames(nm = names(mixed_mod_sp_autocor_list$moran)), function(FD, moran_list) {
  lapply(moran_list[[FD]], function(x) {
    moran_df <- as.data.frame(as.list(x$statistic) )
    moran_df$p_value <- x$p.value
    return(moran_df)
  }) %>% 
    bind_rows(.id = 'predictors')
}, moran_list = mixed_mod_sp_autocor_list$moran) %>% 
  bind_rows(.id = 'outcome') %>% 
  mutate(model = paste(outcome, '~', predictors, sep = " "))



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# summary(mod_rename_list$qEt$`intercept + area + julian_date + temp`)
# 
# 
# summary(mod_rename_list$FRic$`intercept + fish1 + julian_date`)
# summary(mod_rename_list$FDis$`intercept + fish1`)
# summary(mod_rename_list$FDiv$intercept)
# summary(mod_rename_list$FEve$`intercept + fish1 + temp`)
# summary(mod_rename_list$M.prime$`intercept + fish1`)
# summary(mod_rename_list$qEt$`intercept + area + julian_date + temp`)
# summary(mod_rename_list$nsp$`intercept + area + julian_date + temp`)
# summary(mod_rename_list$qDTM$`intercept + fish1`)
# 
# 
# 
# summary(mod_rename_list$FRic$`intercept + fish1`)
# summary(mod_rename_list$FDis$`intercept + fish1`)
# summary(mod_rename_list$FDiv$`intercept + fish1`)
# summary(mod_rename_list$FEve$`intercept + fish1`)
# 
# 
# selection_table_list$FRic
# selection_table_list$FDis
# selection_table_list$FDiv
# selection_table_list$FEve
# 
# summary(mod_rename_list$FRic$`intercept + elev + julian_date`)
# 
# 
# #look at if year shows different relationship
# 
# 
# lapply(alpha_fd_mods$FRic$top_mods)
# 
# mod_var <- paste(rownames(summary(alpha_fd_mods$FRic$top_mods$`11`)$coefficients), collapse = " + ")
# gsub(pattern = '\\(Intercept\\)', 'intercept', mod_var)
# 
# 
# class(alpha_fd_mods$FRic$top_mods)
# 
# 
# selectionTable <- aictab(cand.set = mod_rename_list$FRic)
# 
# library(AER)
# par(mfrow = c(2, 2))
# #,trafo=1
# dispersiontest(mod_rename_list$nsp$`intercept + area`)
#plot(mod_rename_list$nsp$`intercept + area`)

#t.test(x, y, paired = TRUE, alternative = "two.sided")


# ModelAvg <- model.avg(alpha_fd_mods$FDis$dredge_output)
# mA<-summary(ModelAvg) #pulling out model averages
# df1<-as.data.frame(mA$coefmat.full) #selecting full model coefficient averages
# 
# CI <- as.data.frame(confint(ModelAvg, full=T)) # get confidence intervals for full model
# df1$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
# df1$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
# setDT(df1, keep.rownames = "coefficient") #put rownames into column
# names(df1) <- gsub(" ", "", names(df1)) # remove spaces from column headers
# 
# 
# ModelAvg <- model.avg(alpha_fd_mods$FRic$dredge_output)
# mA<-summary(ModelAvg, subset <= 3) #pulling out model averages
# df1<-as.data.frame(mA$coefmat.subset) #selecting full model coefficient averages
# 
# CI <- as.data.frame(confint(ModelAvg)) # get confidence intervals for full model
# df1$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
# df1$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
# #setDT(df1, keep.rownames = "coefficient") #put rownames into column
# #names(df1) <- gsub(" ", "", names(df1)) # remove spaces from column headers
# 
# df1[-1,] %>% 
#   rownames_to_column(var = 'coefficient') %>% 
#   ggplot(aes(x=coefficient, y=Estimate))+ #again, excluding intercept because estimates so much larger
#   geom_hline(yintercept=0, color = "red",linetype="dashed", lwd=1.5)+ #add dashed line at zero
#   geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="blue", #adj SE
#                 width=0, lwd=1.5) +
#   coord_flip()+ # flipping x and y axes
#   geom_point(shape = 21, stroke = 2, size = 5.5) +
#   #geom_point(size=8) + 
#   theme_bw(base_size = 20) + ylab("Coefficient")
# 


# loc_dists <- as.matrix(dist(cbind(fd_ses_df_lake_wider$long, fd_ses_df_lake_wider$lat)))
# loc_dists_inv <- 1/loc_dists
# diag(loc_dists_inv) <- 0
# 
# Moran.I(resid(alpha_fd_mods$FRic$top_mods$`11`), loc_dists_inv)
# 
# ozone.dists.inv <- 1/ozone.dists
# diag(ozone.dists.inv) <- 0
# 
# 
# prac_moran<- moran_wrapper_alt(data = fd_ses_df_lake_wider,
#                       mod = alpha_fd_mods$FDis$top_mods$`5`,
#                       long_col = 'long',
#                       lat_col = 'lat')
# 
# 
# variogram(resid(alpha_fd_mods$FEve$top_mods$`5`) ~ 1, 
#                                      loc = ~long+lat, 
#                                      data = fd_ses_df_lake_wider %>% 
#             mutate(long = (long - mean(long))/sd(long),
#                    lat = (lat - mean(lat))/sd(lat) ))
# 
# hist(as.matrix(dist(fd_ses_df_lake_wider[, c('long', 'lat')])))
# 
# 
# variogram(resid(alpha_fd_mods$FDis$top_mods$`5`) ~ 1, 
#           locations = ~long+lat, 
#           data = fd_ses_df_lake_wider) %>% 
#   ggplot() +
#   geom_point(aes(x = dist, y = gamma))
#   
# 
# alpha_spatial_autocor_list$qDTM
# alpha_spatial_autocor_list$M.prime
# 
# 
# alpha_spatial_autocor_list$FRic$`7`$moran
# 
# par(mfrow = c(2, 2))
# plot(alpha_fd_mods$FDis$top_mods$`5`)
# mtext(paste0("Model 11: ", 'FRic' , ' ~ ', paste((rownames(prac_sum$coefficients)), collapse = ' + ')), side = 3, line = -2, outer = TRUE)
# 
#
# prac_sum <- summary(alpha_fd_mods$FDis$top_mods$`5`)
# 
# paste((rownames(prac_sum$coefficients)), collapse = ', ')
# 
# variogram_plot_list <- lapply(alpha_spatial_autocor_list, function(x) {
#   lapply(x, function(y) {
#     y$variogram %>% 
#       ggplot() +
#       geom_point(aes(x = dist, y = gamma), size = 4, color = 'gray') +
#       ylim(0, max(y$variogram$gamma) ) +
#       xlab('Distance') +
#       ylab('Semivariance') +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#             axis.title = element_text(size = 14))
#       
#   })
# })
# 
# variogram_plot_list$M.prime$`3`
# 
# 
# variogram_plot_list$FRic$`7`
# variogram_plot_list$FDis$`3`
# 
# 
# ggplot(data = alpha_spatial_autocor_list$FRic$`11`$variogram) +
#   geom_point(aes(x = dist, y = gamma), size = 4, color = 'gray') +
#   ylim(0, max(alpha_spatial_autocor_list$FRIc$`11`$variogram$gamma) ) +
#   xlab('Distance') +
#   ylab('Semivariance') +
#   theme_bw() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
#   
# 
# 
# 
# eval_spatial_autocor(data = fd_ses_df_lake_wider, 
#                      mod = alpha_fd_mods$FEve$top_mods$`5`, 
#                      long_col = 'long', 
#                      lat_col = 'lat', 
#                      method = c('moran', 'variogram'))
# 
# 
# variogram(resid(alpha_fd_mods$FEve$top_mods$`5`) ~ 1, 
#           loc = ~long+lat, 
#           data = fd_ses_df_lake_wider %>% 
#             mutate(long = (long - mean(long))/sd(long),
#                    lat = (lat - mean(lat))/sd(lat) )) %>% 
#   ggplot() +
#   geom_point(aes(x = dist, y = gamma))
# 
# 

# 
# # plot_variogram <- function(data,
# #                            point_size = 1,
# #                            plot_title = NULL,
# #                            axis_title_size = 5,
# #                            plot_title_size = 5,
# #                            point_col = 'black') {
# #   
# #   init_plot <- ggplot(data = data) +
# #     geom_point(aes(x = dist, y = gamma), size = point_size, color = point_col) +
# #     ylim(0, max(data$gamma) ) +
# #     xlab('Distance') +
# #     ylab('Semivariance') +
# #     theme_simple() +
# #     theme(axis.title = element_text(size = axis_title_size))
# #   
# #   if (!is.null(plot_title)) {
# #     init_plot <- init_plot +
# #       ggtitle(plot_title) +
# #       theme(plot.title = element_text(size = aplot_title_size))
# #   }
# #   
# #   return(init_plot)
# # }
# 
# 
# 
# 
# 
# export_mod_diagnostics <- function(mod_list,
#                                    mod_set_name,
#                                    export_path,
#                                    file_name,
#                                    diagnostics = c('QQ', 'resid_vs_predict', 'dispersion'),
#                                    panel_width = 6,
#                                    plot_height = 6,
#                                    show_progress = TRUE) {
#   
#   diagnostics <- match.arg(diagnostics, several.ok = TRUE)
#   col_dim <- length(diagnostics)
#   
#   pdf(file = paste0(export_path, file_name, '.pdf'), 
#       width = panel_width, height = plot_height)
#   
#   for (DATASET in names(mod_list)) {
#     
#     for (EVOL_ORIG in names(mod_list[[DATASET]]) ) {
#       
#       for (FD_CONTR in names(mod_list[[DATASET]][[EVOL_ORIG]]) ) {
#         
#         sim_resid <- simulateResiduals(mod_list[[DATASET]][[EVOL_ORIG]][[FD_CONTR]]$mod)
#         
#         par(mfrow = c(1, col_dim), oma = c(1, 1, 2.25, 1))
#         if ('QQ' %in% diagnostics) plotQQunif(sim_resid)
#         if ('resid_vs_predict' %in% diagnostics) plotResiduals(sim_resid, rank = TRUE, quantreg = FALSE, cex = 0.5)
#         if ('dispersion' %in% diagnostics) testDispersion(sim_resid)
#         mtext(paste0(mod_set_name, '; dataset: ', DATASET, '; ', FD_CONTR, ' ~ ',  EVOL_ORIG), 
#               outer = TRUE, side = 3, cex = 1, font = 2)
#         
#         if (show_progress) message(paste0(mod_set_name, '; dataset: ', DATASET, '; ', FD_CONTR, ' ~ ',  EVOL_ORIG))
#       }
#     }
#   }
#   
#   dev.off()
# }
# 
# 
# 
# library(gglm)
# #https://graysonwhite.com/gglm/
# gglm(alpha_fd_mods$FRic$top_mods$`11`)
# gglm(alpha_fd_mods$FEve$top_mods$`5`)
# gglm(alpha_fd_mods$FDiv$top_mods$`9`)
# 
# 
# summary(alpha_fd_mods$FRic$top_mods$`11`)
# summary(alpha_fd_mods$FDis$top_mods$`5`)
# summary(alpha_fd_mods$FEve$top_mods$`5`)
# summary(alpha_fd_mods$FDiv$top_mods$`9`)
# 
# 
# 
# lm(as.formula('fric_ses ~ elev + area + julian_date + time_convert + fish'),
#    data = fd_ses_df_lake_wider, na.action = "na.fail")
# 
# 
# 
# #https://stackoverflow.com/questions/25281739/dredge-function-error-r-package-mumln
# 
# global_mod <- lm(fric_ses ~ elev + area + julian_date + time_convert + fish, 
#                  data = fd_ses_df_lake_wider, na.action = "na.fail")
# 
# fric_dredge <- dredge(global.model = global_mod, beta = c("none"), evaluate = TRUE,
#        rank = "AICc")
# 
# 
# summary(get.models(fric_dredge, subset = delta <= 3)[[1]])
# 
# 
# 
# AICc(lm(fric_ses ~ elev + julian_date + fish, 
#    data = fd_ses_df_lake_wider, na.action = "na.fail"))
# 
# AICc(lm(fric_ses ~ elev + julian_date, 
#         data = fd_ses_df_lake_wider, na.action = "na.fail"))
# 
# AICc(lm(fric_ses ~ fish, 
#         data = fd_ses_df_lake_wider, na.action = "na.fail"))
# 
# AICc(lm(fric_ses ~ elev + julian_date, 
#    data = fd_ses_df_lake_wider, na.action = "na.fail"))
# 
# 
# 
# summary(lm(fric_ses ~ fish, data = fd_ses_df_lake_wider))
# summary(lm(fric_ses ~ elev + julian_date, data = fd_ses_df_lake_wider))
# 
# #https://www.zoology.ubc.ca/~bio501/R/workshops/modelselection.html
# 
# 
# par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
# plot(lm(fric_ses ~ elev + julian_date, data = fd_ses_df_lake_wider))
# par(mfrow=c(1,1))
# 
# 
# global_mod <- lm(fdis_ses ~ elev + area + julian_date + time_convert + fish, 
#                  data = fd_ses_df_lake_wider, na.action = "na.fail")
# 
# fdis_dredge <- dredge(global.model = global_mod, beta = c("none"), evaluate = TRUE,
#                       rank = "AICc")
# 
# 
# get.models(fdis_dredge, 1)

# 
# 
# #https://statmodeling.stat.columbia.edu/2009/07/11/when_to_standar/
# #https://stats.stackexchange.com/questions/29781/when-conducting-multiple-regression-when-should-you-center-your-predictor-varia
# 
#   
# 
# 
# 
# plot(fdis_ses ~ fish, data = fd_ses_df_lake_wider)
# 
# 
# 
# summary(lm(fdis_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# 
# AIC(lm(fdis_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# AIC(lm(fdis_ses ~ fish + julian_date + time_convert, data = fd_ses_df_lake_wider))
# AIC(lm(fdis_ses ~ fish + julian_date, data = fd_ses_df_lake_wider))
# AIC(lm(fdis_ses ~ fish, data = fd_ses_df_lake_wider))
# AIC(lm(fdis_ses ~ fish + area + fish*area, data = fd_ses_df_lake_wider))
# 
# 
# AIC(lm(fric_ses ~ fish + julian_date + time_convert + area + fish*area, data = fd_ses_df_lake_wider))
# AIC(lm(fric_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# AIC(lm(fric_ses ~ fish + julian_date + time_convert, data = fd_ses_df_lake_wider))
# AIC(lm(fric_ses ~ fish + julian_date, data = fd_ses_df_lake_wider))
# AIC(lm(fric_ses ~ fish, data = fd_ses_df_lake_wider))
# AIC(lm(fric_ses ~ fish + area + fish*area, data = fd_ses_df_lake_wider))
# 
# AIC(lm(fdiv_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# AIC(lm(fdiv_ses ~ fish + julian_date + time_convert, data = fd_ses_df_lake_wider))
# AIC(lm(fdiv_ses ~ fish + julian_date, data = fd_ses_df_lake_wider))
# AIC(lm(fdiv_ses ~ fish, data = fd_ses_df_lake_wider))
# 
# AIC(lm(feve_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# AIC(lm(feve_ses ~ fish + julian_date + time_convert, data = fd_ses_df_lake_wider))
# AIC(lm(feve_ses ~ fish + julian_date, data = fd_ses_df_lake_wider))
# AIC(lm(feve_ses ~ fish, data = fd_ses_df_lake_wider))
# 
# summary(lm(fdis_ses ~ fish + julian_date + time_convert + area, data = fd_ses_df_lake_wider))
# 
# summary(lm(fdis_ses ~ fish, data = fd_ses_df_lake_wider))
# plot(lm(fdis_ses ~ fish, data = fd_ses_df_lake_wider))
# 
# fd_ses_df_lake %>% 
#   ggplot() +
#   geom_histogram(aes(x = value), bins = 15) +
#   facet_wrap(~fd_metric)
# 
# fd_ses_df_lake %>% 
#   ggplot() +
#   geom_point(aes(x = value), bins = 15) +
#   facet_wrap(~fd_metric)
# 
# p + facet_grid(vs + carb ~ .)
# 
# 
# chart.Correlation(fd_ses_df_lake_wider[c('temp', 'area', 'julian_date', 'time_convert', 'elev')], 
#                   histogram=TRUE, pch=19,
#                   method = 'spearman')
# 
# #temp, area, date, time, 
# #elev
# 
# 
# cor.test(fd_ses_df_lake_list$fdis_ses$area, 
#          fd_ses_df_lake_list$fdis_ses$temp,
#          method = 'spearman', exact = FALSE)
# 
# 
# fd_ses_df_lake

