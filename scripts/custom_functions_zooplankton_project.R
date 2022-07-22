#########################################
#########################################
### FUNCTIONS FOR ZOOPLANKTON PROJECT ###
#########################################
#########################################

###################################
### FUNCTIONAL SPACE EVALUATION ###
###################################

deviation_index <- function(pcoa,
                            orig_dist,
                            option = c('msd', 'mad'),
                            number_dims = NULL) {
  
  option <- match.arg(option, several.ok = FALSE)
  if (is.null(number_dims)) number_dims <- ncol(pcoa)
  if (number_dims > ncol(pcoa)) stop('number_dims needs to be <= number of columns in pcoa.', call. = FALSE)
  
  max_val <- max(orig_dist) #used in standardization of the PCoA distance matrix
  S <- nrow(pcoa)
  
  #sapply(setNames(seq_len(number_dims), paste0("m_", seq_len(number_dims), "D")), function(dim, pcoa, orig_dist, max_val_orig_distmat, S) {
  results_list <- lapply(stats::setNames(nm = seq_len(number_dims)), function(dim, pcoa, orig_dist, max_val_orig_distmat, S) {
    #5. calculate Euclidean distance matrix on the PCoA
    dist_pcoa <- stats::dist(pcoa[, seq_len(dim), drop = FALSE], method = 'euclidean')
    
    #6. standardize the distance matrix calculated from the PCoA by the maximum value of the original distance matrix
    y_std <- (dist_pcoa/max(dist_pcoa)) * max_val_orig_distmat
    
    #7. calculate MSD/MAD from the original distance matrix and standardized distance matrix calculated from the PCoA
    #return(round( ( (sum((orig_dist - y_std)^2)) / (S*(S - 1)/2) ), 6))
    deviance_transform <- switch(option, mad = {abs(orig_dist - y_std)}, msd = {(orig_dist - y_std)^2})
    #return(round( (sum(deviance_transform) / (S*(S - 1)/2) ), 6))
    
    return(data.frame(dim, round( (sum(deviance_transform) / (S*(S - 1)/2) ), 6)))
    
  }, pcoa = pcoa, orig_dist = orig_dist, max_val_orig_distmat = max_val, S = S)
  
  #process results list
  return(stats::setNames(do.call(rbind, results_list), c('dim', option)))
}



id_elbow <- function(df,
                     xcol, ycol,
                     xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                     round_val = NULL) {
  
  ### 1. CHECKING INPUTS ###
  if (!all(c(xcol, ycol) %in% colnames(df)))
    stop('The xcol and ycol are not in df.', call. = FALSE)
  
  min_max_list <- list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  
  if (any(!sapply(min_max_list, function(x) (is.null(x)) | (is.numeric(x) & length(x) == 1) )))
    stop('xmin, xmax, ymin, and ymax need to be NULL or a numeric of length 1.', call. = FALSE)
  
  if (!(all(sapply(min_max_list[c('xmin', 'ymin')], is.null)) | all(!sapply(min_max_list[c('xmin', 'ymin')], is.null))))
    stop('xmin and ymin either both have to be NULL or both numeric.', call. = FALSE)
  
  if (is.null(xmax) && !is.null(ymax))
    stop("xmax can only be NULL if ymax is also NULL.", call. = FALSE)
  
  if (!is.null(xmax) && is.null(ymin))
    stop("ymax can only be NULL if xmax is also NULL.", call. = FALSE)
  
  if (!(is.null(round_val) | ( is.numeric(round_val) & length(round_val) == 1 )))
    stop('round_val needs to be NULL or a numeric of length 1.', call. = FALSE)
  
  
  ### 2. PROCESSING DATA ###
  df_processed <- df[, c(xcol, ycol)]
  
  if (!is.null(xmin)) df_processed <- rbind(c(xmin, ymin), df_processed)
  if (!is.null(xmax)) df_processed <- rbind(df_processed, c(xmax, ymax))
  
  if (!identical(df_processed[, xcol], order(df_processed[ , xcol]) ))
    df_processed <- df_processed[order(df_processed[ , xcol]), ]
  
  
  ### 3. IDENTIFYING THE ELBOW INFLECTION POINT ###
  
  #3a. use lm to create line between the first and last points
  extreme_points_df <- df_processed[c(1, nrow(df_processed)), ]
  mod <- stats::lm(extreme_points_df[, ycol] ~ extreme_points_df[, xcol])
  
  #3b. constant: predicted y of each x based on the lm performed on the most extreme points
  df_processed$constant <- mod$coef[[1]] + mod$coef[[2]] * df_processed[, xcol]
  
  #3c. benefits: difference between the observed y and the predicted y (constant)
  df_processed$benefits <- df_processed[ , ycol] - df_processed$constant
  
  #3d. id inflectin point --> if the x value is below the x value with maximum benefit, then it is selected
  df_processed$selected <- ifelse(df_processed[, xcol] <= df_processed[which.max(df_processed$benefits), xcol], 'yes', 'no')
  
  
  ### 4. PROCESSING OUTPUT ###
  
  #4a. identify the first and last rows to be included in processed reuslts
  #remove the min and max rows if they were provided in the min/max arguements
  starting_row <- ifelse(is.null(xmin), 1, 2)
  ending_row <- ifelse(is.null(xmax), nrow(df_processed), nrow(df_processed) - 1)
  
  #4b. round the AUC, constant, and benefits columns
  if (!is.null(round_val)) {
    col_eval <- sapply(colnames(df_processed), function(x) x %in% c('AUC', 'constant', 'benefits') )
    df_processed[col_eval] <- lapply(df_processed[col_eval], function(x, round_val) {round(x, round_val)}, round_val = round_val)
  }
  
  #4c. return processed elbow inflection point results
  return(df_processed[starting_row:ending_row, ])
}


coranking_wrapper <- function(orig_dist,
                              input_pcoa,
                              cols = NULL,
                              randomize = FALSE,
                              rand_col_index = NULL) {
  
  if (is.null(cols) ) cols <- ncol(input_pcoa)
  if (ncol(input_pcoa) < cols) cols <- ncol(input_pcoa)
  
  if (isTRUE(randomize)) {
    #set.seed(seed)
    d_dimen <- stats::dist(data.frame(input_pcoa[, seq_len(cols), drop = FALSE][, -rand_col_index, drop = FALSE], rand_axis = sample(input_pcoa[, rand_col_index])))
  } else {
    d_dimen <- stats::dist(input_pcoa[, seq_len(cols), drop = FALSE])
  }
  
  co_rank <- suppressWarnings(coRanking::coranking(orig_dist, d_dimen, input_Xi = "dist"))
  nx <- coRanking::R_NX(co_rank)
  
  return(coRanking::AUC_ln_K(nx))
  #return(list(corank = co_rank,
  #            nx = nx,
  #            auc = coRanking::AUC_ln_K(nx) ) )
  
}



corank_dim_eval <- function(pcoa,
                            orig_dist,
                            number_dims = NULL,
                            rand_rep = 99) {
  
  if (is.null(number_dims)) number_dims <- ncol(pcoa)
  if (number_dims > ncol(pcoa)) stop('number_dims needs to be <= number of columns in pcoa', call. = FALSE)
  
  compute_auc_and_nullmod <- lapply(seq_len(number_dims), function(dim, pcoa, orig_dist, rand_rep) {
    
    dat <- data.frame(dim_range = dim,
                      auc = coranking_wrapper(orig_dist = orig_dist, input_pcoa = pcoa, cols = dim, randomize = FALSE))
    
    if (dim == 1) {
      rand_vec <- NA
    } else {
      rand_vec <- replicate(rand_rep, coranking_wrapper(orig_dist = orig_dist,
                                                        input_pcoa = pcoa,
                                                        cols = dim,
                                                        randomize = TRUE,
                                                        rand_col_index = dim),
                            simplify = 'vector')
    }
    
    return(list(obs = dat, rand = rand_vec))
    
  }, pcoa = pcoa, orig_dist = orig_dist, rand_rep = rand_rep)
  
  auc_ses_pval <- do.call(rbind, lapply(seq_len(number_dims), function(dim, compute_auc_and_nullmod) {
    
    #AUC_obs <- compute_auc_and_nullmod[[k]]$obs$auc
    
    if (dim == 1) {
      return(data.frame(auc = compute_auc_and_nullmod[[dim]]$obs$auc, ses = NA, pval = NA))
    } else {
      ses <- (compute_auc_and_nullmod[[dim]]$obs$auc - mean(compute_auc_and_nullmod[[dim]]$rand))/stats::sd(compute_auc_and_nullmod[[dim]]$rand)
      # p-value = proportion of null value inferior to obs Beta (+1)
      pval <- sum(compute_auc_and_nullmod[[dim]]$rand > compute_auc_and_nullmod[[dim]]$obs$auc)/(length(compute_auc_and_nullmod[[dim]]$rand) + 1)
    }
    
    return(data.frame(auc = compute_auc_and_nullmod[[dim]]$obs$auc, ses = ses, pval = pval))
    
  }, compute_auc_and_nullmod = compute_auc_and_nullmod))
  
  elbow_eval <- id_elbow(df = data.frame(dims = seq_len(number_dims), auc = auc_ses_pval$auc),
                         xcol = 'dims', ycol = 'auc',
                         xmin = 0, xmax = NULL, ymin = 0, ymax = NULL,
                         round_val = 6)
  
  return(
    data.frame(dim = seq_len(number_dims),
               auc = elbow_eval$auc,
               benefit_auc_elbow = elbow_eval$benefits,
               selected_by_auc_elbow = elbow_eval$selected,
               ses = auc_ses_pval$ses,
               pvalue = auc_ses_pval$pval)
  )
}


eval_func_space <- function(trait_data,
                            dim = 5,
                            dist_method = c("euclidean", "gowdis", "gawdis"),
                            dist_arg_list = NULL,
                            eval_option = c('corank', 'msd', 'mad'),
                            correction_method = c("quasieuclid", "cailliez", "lingoes", "sqrt"),
                            rand_rep = 99,
                            show_messages = TRUE) {
  
  #========================================================
  #== Preparatory steps for functional space evaluation ===
  #========================================================
  
  ### CHECKING AND PROCESSING INPUTS ###
  input_type <- switch(class(trait_data),
                       data.frame = 'trait',
                       dist = 'dist',
                       stop('`trait_data` needs to be a dataframe of traits (class "data.frame") or species distances (class "dist").', call. = FALSE))
  
  #trait_data
  if (nrow(as.matrix(trait_data)) < 2)
    stop('There must be at least 2 species in `trait_data`.', call. = FALSE)
  
  #eval_option
  eval_option <- match.arg(eval_option, several.ok = TRUE)
  eval_option <- sort(unique(eval_option))
  
  #dim
  if (dim <= 0)
    stop('`dim` must be greater than 0.', call. = FALSE)
  
  #arguments to check only if dissimilarity is created internally
  if (input_type == 'trait') {
    dist_method <- match.arg(dist_method, several.ok = FALSE)
    
    if (!(class(dist_arg_list) %in% c('NULL', 'list')))
      stop('`dist_arg_list` needs to be a list of function arguments or NULL. Default is NULL.', call. = FALSE)
    
    #if (!is.logical(std) | length(std) != 1)
    #  stop('"std" needs to be a single logical value (TRUE or FALSE). Default is TRUE.', call. = FALSE)
  }
  
  #show_messages
  #if (!is.logical(show_messages) | length(show_messages) != 1 )
  if (!(isTRUE(show_messages) | isFALSE(show_messages)) )
    stop('`show_messages` needs to be a single logical value (TRUE or FALSE). Default is TRUE.', call. = FALSE)
  
  #correction_method
  correction_method <- match.arg(correction_method, several.ok = TRUE)
  
  #rand_rep
  if (!is.numeric(rand_rep) | length(rand_rep) != 1)
    stop('`rand_rep` accepts only a single, numeric value. Default is 99.', call. = FALSE)
  
  #if std == TRUE and input_type == 'trait' standardize traits in the trait dataset
  #if (isTRUE(std) & input_type == 'trait') {
  #  trait_data[sapply(trait_data, is.numeric)] <- lapply(trait_data[sapply(trait_data, is.numeric)], standardize_func)
  #}
  
  
  ### WORKING WITH THE DISSIMILARITY MATRIX ###
  if (input_type == 'trait') {
    if (isTRUE(show_messages)) message("Calculating dissimilarity matrix.")
    
    #1. calculate distance matrix from the trait data
    dist_mat <-  dist_wrapper(df = trait_data,
                              method = switch(dist_method, euclidean = {'dist'}, gowdis = {'gowdis'}, gawdis = {'gawdis'}),
                              arg_list = dist_arg_list)
    
    distance_matrix_info <- list(dist_mat = list(original_dist = dist_mat))
    #distance_matrix_info <- list(dist_mat = list(original_dist = dist_mat$dist)) #older dist_wrapper function
    
    if (is.null(dist_arg_list)) {
      distance_matrix_info[['dist_args']] <- 'NUll input; function defaults used for dissimilarity matrix construction'
    } else  distance_matrix_info[['dist_args']] <- dist_arg_list
    
  } else { #input_type == dist
    distance_matrix_info <- list(dist_mat = list(original_dist = trait_data))
  }
  
  
  #check to make sure that dissimilarity matrix doesn't contain any NA entries
  if (any(is.na(distance_matrix_info[["dist_mat"]][["original_dist"]])))
    stop("Unable to proceed with functional space evaluation because of missing values in dissimilarity matrix.", call. = FALSE)
  
  ####################
  #return(distance_matrix_info[["dist_mat"]][["original_dist"]])
  ####################
  
  #message about whether the dissimilarity matrix is Euclidean (this is included as an output)
  distance_matrix_info[["euclidean_eval"]] <- ifelse(suppressWarnings(ade4::is.euclid(distance_matrix_info[["dist_mat"]][["original_dist"]])),
                                                     "Original dissimilarity matrix is Euclidean",
                                                     "Original dissimilarity matrix is not Euclidean")
  
  #if the original dissimilarity matrix is not Euclidean, correct the dissimilarity matrix
  if (suppressWarnings(ade4::is.euclid(distance_matrix_info[["dist_mat"]][["original_dist"]]))) {
    transform_vec <- 'untransformed'
    distance_matrix_info[["dist_mat"]][[transform_vec]] <- distance_matrix_info[["dist_mat"]][["original_dist"]]
  } else {
    transform_vec <- sort(unique(correction_method))
    
    for (transformation in transform_vec) {
      distance_matrix_info[["dist_mat"]][[transformation]] <- suppressWarnings(
        switch(transformation,
               quasieuclid = {ade4::quasieuclid(distance_matrix_info[["dist_mat"]][["original_dist"]])},
               cailliez = {ade4::cailliez(distance_matrix_info[["dist_mat"]][["original_dist"]])},
               lingoes = {ade4::lingoes(distance_matrix_info[["dist_mat"]][["original_dist"]])},
               sqrt = {sqrt(distance_matrix_info[["dist_mat"]][["original_dist"]])},
               stop("Enter at least one valid correction option: quasieuclid, cailliez, lingoes, and/or sqrt.", call. = FALSE))
      )
    }
  }
  
  
  ### CREATING LISTS TO STORE RESULTS ###
  PCoA_list <- list()
  func_eval_list <- list()
  
  
  
  #==================================
  #== Functional space evaluation ===
  #==================================
  
  if (isTRUE(show_messages)) message("Finished evaluating: ", appendLF = FALSE)
  
  #id_euclid_vec records whether a PCoA was successfully created from a particular transformation
  id_euclid_vec <- stats::setNames(rep(FALSE, length(transform_vec)), transform_vec)
  
  for (transformation in transform_vec) {
    
    #3. compute PCoA on corrected distance matrix and add to PCoA_list
    #pcoa_temp <- NA
    #suppressWarnings(try(pcoa_temp <- ade4::dudi.pco(distance_matrix_info[["dist_mat"]][[transformation]], scannf = FALSE, full = TRUE), silent = TRUE))
    
    pcoa_temp <- suppressWarnings(tryCatch(ade4::dudi.pco(distance_matrix_info[["dist_mat"]][[transformation]], scannf = FALSE, full = TRUE),
                                           error = function(condition) NA))
    
    if (suppressWarnings(ade4::is.euclid(distance_matrix_info[["dist_mat"]][[transformation]]) & !identical(NA, pcoa_temp) ))  {
      
      #change id_euclid_vec to TRUE if PCoA was successfully created
      id_euclid_vec[transformation] <- TRUE
      
      PCoA_list[[transformation]] <- pcoa_temp$li
      
      #4. determine whether dim or dimensions of PCoA is less (the smaller of these two values will be the number of axes used in functional space evaluation)
      if (dim > ncol(PCoA_list[[transformation]])) {
        warning('`dim` (', dim, ') is greater than the dimensionality of the PCoA (', ncol(PCoA_list[[transformation]]), ') for the ', transformation, ' transformation.\nEvaluating functional space over the available dimensions of the PCoA.', call. = FALSE)
        dim <- ncol(PCoA_list[[transformation]])
      }
      
      
      for (method in eval_option) {
        
        if (method == 'corank') {
          func_eval_list[[method]][[transformation]] <- corank_dim_eval(pcoa = PCoA_list[[transformation]],
                                                                        orig_dist = distance_matrix_info[["dist_mat"]][["original_dist"]],
                                                                        number_dims = dim,
                                                                        rand_rep = rand_rep)
        } else {
          func_eval_list[[method]][[transformation]] <- deviation_index(pcoa = PCoA_list[[transformation]],
                                                                        orig_dist = distance_matrix_info[["dist_mat"]][["original_dist"]],
                                                                        option = method,
                                                                        number_dims = dim)
        }
      }
    }
    
    if (isTRUE(show_messages)) {
      if (transformation == transform_vec[1]) {
        message(transformation, appendLF = FALSE)
      } else if (transformation != transform_vec[length(transform_vec)]) {
        message(', ', transformation, appendLF = FALSE)
      } else {
        message(', ', transformation, '.', appendLF = TRUE)
      }
    }
    
  }
  
  
  
  #=============================
  #== Preparation of outputs ===
  #=============================
  
  if (all(!id_euclid_vec)) {
    func_eval_list <- 'No functional space evaluation conducted.'
  } else {
    func_eval_list <- lapply(func_eval_list, function(list_element) {
      list1 <- listname2column(list_element, column_name = 'correction')
      return(do.call(rbind, unname(list1)))
    } )
    
    
    # for (method in eval_option[eval_option != 'corank']) {
    #   euclid_deviation <- func_eval_list[[method]][id_euclid_vec]
    #
    #   #number of PCoA dimensions that MSD/MAD was calculated over for each correction and/or distance matrix method
    #   k_dim <- sapply(euclid_deviation, length)
    #
    #   #if the number of dimensions that MSD/MAD was calculated over is not uniform, add NAs to MSD/MAD calculation vectors that
    #   #have less than the max calculations so that all calculation vectors have the same length. This is necessary for
    #   #correctly combining the calculations into a single dataframe.
    #   if (length(unique(k_dim) ) != 1) {
    #     euclid_deviation <- lapply(euclid_deviation, function(x, max_k) {c(x, rep(NA, max_k - length(x)))}, max_k = max(k_dim))
    #   }
    #
    #   #collapse the list of MSD/MAD calculations to a dataframe
    #   deviation_dataframe <- data.frame(do.call(cbind, euclid_deviation))
    #   func_eval_list[[method]] <- cbind(dim = seq_len(nrow(deviation_dataframe)), deviation_dataframe)
    #
    #   #func_eval_list[[method]][['results']] <- cbind(dim = seq_len(nrow(deviation_dataframe)), deviation_dataframe)
    #   #func_eval_list[[method]] <- func_eval_list[[method]][!(names(func_eval_list[[method]]) %in% names(id_euclid_vec)[id_euclid_vec])]
    # }
  }
  
  return(
    list(func_eval = func_eval_list,
         supplementary_info = list(PCoA_list = PCoA_list,
                                   distance_matrix_info = distance_matrix_info,
                                   euclid_info = id_euclid_vec))
  )
}


parse_deviation_output <- function(deviation_df,
                                   max_dim = NULL,
                                   deviation_name) {
  
  if (!(deviation_name %in% colnames(deviation_df)))
    stop('The name of the deviationn column you provided does not exist in deviation_df', call. = FALSE)
  
  if (!is.null(max_dim)) {
    
    best_dim_df <- 'max_dim_results'
    
    if (max_dim > max(deviation_df$dim)) {
      warning('`max_dim` is greater than the dimensionality of `deviation_df`.\nAll dimensions will be considered when identifying optimal dimension/correction combination.', call. = FALSE)
      max_dim <- max(deviation_df$dim)
    }
  } else {
    best_dim_df <- 'all_dim_results'
  }
  
  results_list <- list()
  
  #all dim results
  results_list[['all_dim_results']] <- deviation_df %>%
    dplyr::group_by(.data$correction) %>%
    dplyr::slice_min(.data[[deviation_name]], n = 1, with_ties = FALSE) %>%
    as.data.frame()
  
  if (!is.null(max_dim)) {
    results_list[['max_dim_results']] <- deviation_df %>%
      dplyr::filter(.data$dim <= max_dim) %>%
      dplyr::group_by(.data$correction) %>%
      dplyr::slice_min(.data[[deviation_name]], n = 1, with_ties = FALSE) %>%
      as.data.frame()
  }
  
  results_list[['best_dim_info']] <- results_list[[best_dim_df]] %>%
    dplyr::slice_min(.data[[deviation_name]], n = 1, with_ties = FALSE)
  
  return(results_list)
}


parse_corank_output <- function(corank_df,
                                max_dim = NULL,
                                upper_auc_thresh = 0.7,
                                lower_auc_thresh = NULL) {
  
  #if (!is.null(ignore)) corank_list <- corank_list[!(names(corank_list) %in% ignore)]
  #corank_list <- corank_list[sapply(corank_list, is.data.frame)]
  
  corank_list <- split(corank_df, corank_df$correction)
  
  if (is.null(max_dim)) {
    #if max dim is NULL, use all possible dimensions
    max_dim_vec <- sapply(corank_list, nrow)
  } else if (any(max_dim > sapply(corank_list, nrow))) {
    warning('`max_dim` is greater than the dimensionality of least one of the corank results.\nUsing the maximum dimensionality for the corrections with dimensionality less than `max_dim`.', call. = FALSE)
    max_dim_vec <- sapply(corank_list, nrow)
    max_dim_vec[max_dim_vec > max_dim] <- max_dim
  } else {
    #use max_dim for all corrections if max_dim is less than the dimensionality of all corrections
    max_dim_vec <- stats::setNames(rep(max_dim, length(corank_list)), names(corank_list))
  }
  
  #1. identify the highest dim that the elbow method identified
  #elbow_selection <- lapply(corank_list, function(df) df[max(which(df$selected_by_auc_elbow == 'yes')),] )
  elbow_selection <- lapply(corank_list, function(df) df[df$dim == max(df[df$selected_by_auc_elbow == 'yes',]$dim),])
  
  #2. check if elbow selection is below the max_dim and identify the best dim <= max_dim if not
  elbow_selection_updated <- lapply(stats::setNames(nm = names(elbow_selection)),
                                    function(correction, elbow_selection, corank_list, max_dim_vec) {
                                      
                                      if (elbow_selection[[correction]]$dim <= max_dim_vec[correction]) {
                                        return(elbow_selection[[correction]])
                                      } else {
                                        subset_df <- corank_list[[correction]][corank_list[[correction]]$dim %in% seq_len(max_dim_vec[correction]),]
                                        return(subset_df[subset_df$auc == max(subset_df$auc),][1,])
                                      }
                                    }, corank_list = corank_list, elbow_selection = elbow_selection, max_dim_vec = max_dim_vec)
  
  #3. if the updated elbow selection is below the max dims and AUC < upper_auc_thresh, increase the number of dims until you reach the max or you each upper_auc_thresh
  final_selection <- lapply(stats::setNames(nm = names(elbow_selection_updated)),
                            function(correction, elbow_selection, corank_list, max_dim_vec) {
                              if (elbow_selection[[correction]]$dim < max_dim_vec[correction]) {
                                subset_df <- corank_list[[correction]][corank_list[[correction]]$dim %in% elbow_selection[[correction]]$dim:max_dim_vec[correction],]
                                return(subset_df[subset_df$auc == max(subset_df$auc),][1,])
                              } else {
                                return(elbow_selection[[correction]])
                              }
                            }, corank_list = corank_list, elbow_selection = elbow_selection_updated, max_dim_vec = max_dim_vec)
  
  #final_selection_df <- data.frame(correction = names(final_selection), do.call(rbind, unname(final_selection) ))
  final_selection_df <- do.call(rbind, unname(final_selection) )
  if (!is.null(lower_auc_thresh)) final_selection_df$meets_lower_auc_thresh <- final_selection_df$auc >= lower_auc_thresh
  
  return(list(elbow_selection = do.call(rbind, unname(elbow_selection)),
              processed_results = final_selection_df,
              best_dim_info = final_selection_df[final_selection_df$auc == max(final_selection_df$auc),])
  )
  
}


# parse_deviation_output <- function(deviation_df,
#                                    max_dim = 6,
#                                    deviation_name = 'deviation') {
#
#   #if (!is.null(ignore)) deviation_df <- deviation_df[, !(colnames(deviation_df) %in% ignore)]
#
#   if (max_dim > nrow(deviation_df))
#     warning('`max_dim` is greater than the dimensionality (number of rows) of `deviation_df`.\nAll dimensions will be considered when identifying optimal dimension/correction combination.', call. = FALSE)
#
#   deviation_df <- deviation_df[, sapply(deviation_df, function(x) !all(is.na(x)))]
#
#   subset_df <- deviation_df[seq_len(min(max_dim, nrow(deviation_df))),]
#
#   results_list <- lapply(stats::setNames(seq_len(2), c('all_dim_results', 'max_dim_results')), function(x, deviation_name) {
#     results_df <- data.frame(colnames(deviation_df)[colnames(deviation_df) != 'dim'], NA, NA)
#     colnames(results_df) <- c('correction', 'dim', deviation_name)
#     return(results_df)
#   }, deviation_name = deviation_name)
#
#
#   for (correct in colnames(deviation_df)[colnames(deviation_df) != 'dim']) {
#
#     #full results
#     low_dev_ind <- which(deviation_df[, correct] == min(deviation_df[, correct], na.rm = TRUE)  )[1]
#
#     results_list[['all_dim_results']][which(results_list[['all_dim_results']]$correction == correct), 'dim'] <- deviation_df[low_dev_ind,]$dim
#     results_list[['all_dim_results']][which(results_list[['all_dim_results']]$correction == correct), deviation_name] <- deviation_df[low_dev_ind, correct]
#
#     #subset results
#     low_dev_ind <- which(subset_df[, correct] == min(subset_df[, correct], na.rm = TRUE))[1]
#     results_list[['max_dim_results']][which(results_list[['max_dim_results']]$correction == correct), 'dim'] <- subset_df[low_dev_ind,]$dim
#     results_list[['max_dim_results']][which(results_list[['max_dim_results']]$correction == correct), deviation_name] <- subset_df[low_dev_ind, correct]
#   }
#
#
#   #best overall dimension for each correction
#   results_list[['best_dim_info']] <- results_list[['max_dim_results']][results_list[['max_dim_results']][, deviation_name] == min(results_list[['max_dim_results']][, deviation_name]),]
#
#   return(results_list)
#
# }

# parse_corank_output <- function(corank_list,
#                                 max_dim = 6,
#                                 upper_auc_thresh = 0.7,
#                                 lower_auc_thresh = NULL) {
#
#   #if (!is.null(ignore)) corank_list <- corank_list[!(names(corank_list) %in% ignore)]
#   corank_list <- corank_list[sapply(corank_list, is.data.frame)]
#
#   if (any(max_dim > sapply(corank_list, nrow)))
#     stop('`max_dim` is greater than the dimensionality of least one of the corank results', call. = FALSE)
#
#   #1. identify the highest dim that the elbow method identified
#   elbow_selection <- lapply(corank_list, function(df) df[max(which(df$selected_by_auc_elbow == 'yes')),] )
#
#   #2. check if elbow selection is below the max_dim and identify the best dim <= max_dim if not
#   elbow_selection_updated <- lapply(stats::setNames(nm = names(elbow_selection)),
#                                     function(correction, elbow_selection, corank_list, max_dim) {
#
#                                       if (elbow_selection[[correction]]$dim <= max_dim) {
#                                         return(elbow_selection[[correction]])
#                                       } else {
#                                         subset_df <- corank_list[[correction]][corank_list[[correction]]$dim %in% seq_len(max_dim),]
#                                         return(subset_df[subset_df$auc == max(subset_df$auc),][1,])
#                                       }
#                                     }, corank_list = corank_list, elbow_selection = elbow_selection, max_dim = max_dim)
#
#   #3. if the updated elbow selection is below the max dims and AUC < upper_auc_thresh, increase the number of dims until you reach the max or you each upper_auc_thresh
#   final_selection <- lapply(stats::setNames(nm = names(elbow_selection_updated)),
#                             function(correction, elbow_selection, corank_list, max_dim) {
#                               if (elbow_selection[[correction]]$dim < max_dim) {
#                                 subset_df <- corank_list[[correction]][corank_list[[correction]]$dim %in% elbow_selection[[correction]]$dim:max_dim,]
#                                 return(subset_df[subset_df$auc == max(subset_df$auc),][1,])
#                               } else {
#                                 return(elbow_selection[[correction]])
#                               }
#                             }, corank_list = corank_list, elbow_selection = elbow_selection_updated, max_dim = max_dim)
#
#   final_selection_df <- data.frame(correction = names(final_selection), do.call(rbind, final_selection))
#   if (!is.null(lower_auc_thresh)) final_selection_df$meets_lower_auc_thresh <- final_selection_df$auc >= lower_auc_thresh
#
#   return(list(elbow_selection = data.frame(correction = names(elbow_selection), do.call(rbind, elbow_selection)),
#               processed_results = final_selection_df,
#               best_dim_info = final_selection_df[final_selection_df$auc == max(final_selection_df$auc),])
#   )
# }

#' Parse functional space quality information
#'
#' Parses the output of \code{\link{eval_func_space}} to determine which combination of PCoA dimensionality and
#' Euclidean correction method (if necessary) that forms the optimal functional space.
#'
#' @param eval_func_space_output The functional space evaluation results from \code{\link{eval_func_space}} to be
#'    parsed.
#' @param max_dim The maximum number of dimensions considered when identifying the optimal functional
#'    space. This argument can be useful if future methods performed on the functional space can only
#'    handle a particular number of dimensions, and thus you would want to identify the optimal functional
#'    space within the dimensionality constraint.
#' @param lower_auc_thresh (optional) A lower threshold for the area under the curve (AUC) calculated as part of the
#'    co-ranking matrix method. A logical column is added to the output where functional spaces with AUC
#'    greater or equal to \code{lower_auc_thresh} have a value of \code{TRUE} and functional spaces with AUC below
#'    \code{lower_auc_thresh} have a value of \code{FALSE}. This argument is included because it is important to
#'    consider not only the relative performance of functional spaces (e.g. which functional space is the best
#'    compared to others), but also their absolute performance (e.g. is the best functional space actually high
#'    quality?). In this context, \code{lower_auc_thresh} can be set specify a lower quality threshold below which
#'    functional spaces have unacceptably low quality. For example,
#'    \href{https://doi.org/10.1111/ele.13778}{Mouillot et al. (2021)} suggest that an AUC value below 0.5 indicates
#'    problematically low quality. If you do not want consider a lower AUC threshold, set \code{lower_auc_thresh} to
#'    \code{NULL} (the default).
#' @param upper_auc_thresh A threshold used to identify high quality functional spaces.
#'
#' @return returns some stuff
#' @export
#'
#' @examples sum(1)
parse_eval_func_space <- function(eval_func_space_output,
                                  max_dim,
                                  lower_auc_thresh = NULL,
                                  upper_auc_thresh = 0.7) {
  
  ### Checking inputs ###
  arg_list <- stats::setNames(list(max_dim, upper_auc_thresh),
                              c('max_dim', 'upper_auc_thresh'))
  for (argument in c('max_dim', 'upper_auc_thresh')) {
    if (!is.numeric(arg_list[[argument]]) | length(arg_list[[argument]]) != 1)
      stop(argument, ' needs to be a single, numeric value.', call. = FALSE)
  }
  
  if ( ( !is.numeric(lower_auc_thresh) | length(lower_auc_thresh) != 1 ) & !is.null(lower_auc_thresh))
    stop('`lower_auc_thresh` needs to be a single numeric value or NULL.', call. = FALSE)
  
  if (!is.null(lower_auc_thresh)) {
    if (lower_auc_thresh < 0 | lower_auc_thresh > 1) {
      warning('AUC is constrained between 0 and 1. Changing `lower_auc_thresh` (currently ', lower_auc_thresh, ') to 0.', call. = FALSE)
      lower_auc_thresh <- 0
    }
    
    if (lower_auc_thresh < 0 | lower_auc_thresh > 1) {
      warning('AUC is constrained between 0 and 1. Changing `upper_auc_thresh` (currently ', upper_auc_thresh, ') to 1.', call. = FALSE)
      upper_auc_thresh <- 1
    }
    
    if (upper_auc_thresh < lower_auc_thresh) {
      warning('`lower_auc_thresh` (', lower_auc_thresh, ') is bigger than `upper_auc_thresh` (', upper_auc_thresh, '). Changing `lower_auc_thresh` to equal `upper_auc_thresh`.', call. = FALSE)
      lower_auc_thresh <- upper_auc_thresh
    }
  }
  
  
  ### Processing functional space evaluation results ###
  
  parsed_func_space_eval <- lapply(stats::setNames(nm = names(eval_func_space_output)), function(eval_method, results_list) {
    
    switch(eval_method,
           corank = parse_corank_output(corank_df = results_list[['corank']],
                                        max_dim = max_dim,
                                        upper_auc_thresh = upper_auc_thresh,
                                        lower_auc_thresh = lower_auc_thresh),
           mad = parse_deviation_output(deviation_df = results_list[[eval_method]],
                                        max_dim = max_dim,
                                        deviation_name = eval_method),
           msd = parse_deviation_output(deviation_df = results_list[[eval_method]],
                                        max_dim = max_dim,
                                        deviation_name = eval_method),
           stop('The elements of `eval_func_space_output` should be named by the functional evaluation method (msd, mad, and/or corank).', call. = FALSE))
    
  }, results_list = eval_func_space_output)
  
  return(list(parsed_results = parsed_func_space_eval,
              parse_criteria = stats::setNames(list(max_dim, lower_auc_thresh, upper_auc_thresh),
                                               c('max_dim', 'lower_auc_thresh', 'upper_auc_thresh')))
  )
}


plot_func_space_eval <- function(eval_func_space_output,
                                 max_dim,
                                 upper_auc_thresh = 0.7,
                                 label_top_spaces = TRUE,
                                 methods_color = NULL,
                                 correction_color = NULL,
                                 text_base_size = 12) {
  
  ### CHECK INPUTS ###
  if (!(isTRUE(label_top_spaces) | isFALSE(label_top_spaces)) )
    stop('`label_top_spaces` needs to be a single logical value (TRUE or FALSE). Default is TRUE.', call. = FALSE)
  
  
  ### PROCESS DATASETS ###
  # initial_process <- lapply(stats::setNames(nm = names(eval_func_space_output)), function(eval_method, results_list) {
  #
  #   switch(eval_method,
  #          corank = {
  #            # processed_corank_list <- lapply(stats::setNames(nm = names(results_list$corank)), function(x, corank_list) {
  #            #   corank_list[[x]]$Correction <- x
  #            #   return(corank_list[[x]])
  #            # }, corank_list = results_list$corank)
  #
  #            # dplyr::bind_rows(processed_corank_list) %>%
  #            #   dplyr::select(c('dim', 'auc', 'Correction') ) %>%
  #            #   dplyr::mutate(Correction = factor(.data$Correction, levels = sort(unique(.data$Correction)) ))
  #
  #            results_list$corank %>%
  #              dplyr::select(c('dim', 'auc', 'Correction') ) %>%
  #              dplyr::mutate(Correction = factor(.data$Correction, levels = sort(unique(.data$Correction)) ))
  #
  #          },
  #          mad = {
  #            results_list[[eval_method]] %>%
  #              tidyr::pivot_longer(!.data$dim, names_to = "Correction", values_to = 'mad') %>%
  #              dplyr::mutate(Correction = factor(.data$Correction, levels = sort(unique(.data$Correction)))) %>%
  #              as.data.frame()
  #
  #          },
  #          msd = {
  #            results_list[[eval_method]] %>%
  #              tidyr::pivot_longer(!.data$dim, names_to = "Correction", values_to = 'msd') %>%
  #              dplyr::mutate(Correction = factor(.data$Correction, levels = sort(unique(.data$Correction)))) %>%
  #              as.data.frame()
  #          },
  #          stop('The elements of `eval_func_space_output` should be named by the functional evaluation method (msd, mad, and/or corank).', call. = FALSE))
  #
  # }, results_list = eval_func_space_output)
  
  eval_func_space_output_vis <- lapply(eval_func_space_output, function(x) x %>% dplyr::rename('Correction' = 'correction'))
  
  #initiate the plot(s)
  #plot_list <- lapply(initial_process, function(x) x %>% ggplot2::ggplot())
  plot_list <- lapply(eval_func_space_output_vis, function(x) x %>% ggplot2::ggplot())
  
  if (isTRUE(label_top_spaces)) {
    #pfse_vis <- parse_eval_func_space(eval_func_space_output = eval_func_space_output,
    pfse_vis <- parse_eval_func_space(eval_func_space_output = eval_func_space_output,
                                      max_dim = max_dim,
                                      lower_auc_thresh = NULL,
                                      upper_auc_thresh = upper_auc_thresh)
    
    
    #pfse_vis_extract <- lapply(pfse_vis$parsed_results, function(x) x$best_dim_info)
    pfse_vis_extract <- do.call(rbind, lapply(stats::setNames(nm = names(pfse_vis$parsed_results)), function(x, results) {
      
      #method_name <- switch(x, corank = 'co-rank', toupper(x))
      updated_method <- switch(x, corank = 'co-rank', 'MSD/MAD')
      
      results[[x]]$best_dim_info %>%
        dplyr::rename('Correction' = 'correction') %>%
        dplyr::select(c('Correction', 'dim')) %>%
        dplyr::mutate(old_method = x,
                      Method = updated_method)
    }, results = pfse_vis$parsed_results))
    
    pfse_vis_extract2 <- pfse_vis_extract[!duplicated(pfse_vis_extract$Method),]
    
    plot_list <- lapply(stats::setNames(nm = names(plot_list)), function(method, plot_list, results, eval_func_space_output_vis) {
      y_var <- switch(method, corank = 'auc', method)
      
      plot_list[[method]] +
        ggplot2::geom_point(data = dplyr::left_join(results, eval_func_space_output_vis[[method]],
                                                    by = c('Correction', 'dim')),
                            ggplot2::aes(x = .data[['dim']], y = .data[[y_var]], color = .data[['Method']]),
                            size = 5.5, pch = 21, stroke = 2, fill = NA)
    }, plot_list = plot_list, results = pfse_vis_extract2, eval_func_space_output_vis = eval_func_space_output_vis)
  }
  
  plot_list_final <- lapply(stats::setNames(nm = names(plot_list)), function(eval_method, plot_list, processed_data) {
    
    if (eval_method == 'corank') {
      plot_list[[eval_method]] +
        ggplot2::geom_point(ggplot2::aes(x = .data$dim, y = .data[['auc']], fill = .data$Correction),
                            pch = 21, size = 4, color = "#d8e2dc", stroke = 0.1) +
        ggplot2::xlab('PCoA dimensionality') + ggplot2::ylab('AUC') +
        ggplot2::ggtitle('Co-ranking matrix') +
        ggplot2::scale_x_continuous(breaks = seq_len(max(processed_data[[eval_method]]$dim))) +
        ggplot2::ylim(0, 1)
    } else { #msd/mad
      plot_list[[eval_method]] +
        ggplot2::geom_point(ggplot2::aes(x = .data$dim, y = .data[[eval_method]], fill = .data$Correction),
                            pch = 21, size = 4, color = "#d8e2dc", stroke = 0.1) +
        ggplot2::xlab('PCoA dimensionality') + ggplot2::ylab(toupper(eval_method)) +
        ggplot2::ggtitle(toupper(eval_method)) +
        ggplot2::scale_x_continuous(breaks = seq_len(max(processed_data[[eval_method]]$dim))) +
        ggplot2::ylim(0, max(processed_data[[eval_method]][, eval_method]))
    }
  }, plot_list = plot_list, processed_data = eval_func_space_output_vis)
  
  return(plot_list_final)
}



########################################
### CALCULATING FUNCTIONAL DIVERSITY ###
########################################

FDis_PCoA_input <- function(d, tol = 1e-07) {
  n <- attr(d, "Size")
  
  if (any(is.na(d) ) ) stop("NA's in the distance matrix.","\n")
  A <- matrix(0, ncol = n, nrow = n)
  
  A[row(A) > col(A)] <- -0.5 * d^2
  
  A <- A + t(A)
  
  # gower's double-centering
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  
  vectors <- e$vectors
  eig <- e$values
  
  # check if 'd' is Euclidean or not
  w0 <- eig[n] / eig[1]
  
  if (w0 > -tol) {
    r <- sum(eig > (eig[1] * tol))
  } else {
    r <- length(eig)
  }
  
  vectors <- vectors[, seq_len(r), drop = FALSE] %*% diag(sqrt(abs(eig <- eig[seq_len(r)])), r)
  
  rownames(vectors) <- attr(d, "Labels") #colnames(as.matrix(d))
  
  pos <- eig > 0
  
  return(list(vectors = vectors, pos = pos) )
}





fdisp_alt <- function(d, a, input, tol = 1e-07) {
  
  
  if (input == "distance") {
    
    if(!inherits(d, "dist"))
      stop("'d' must be a 'dist' object.")
    n <- attr(d, "Size")
    
    if (is.null(attr(d, "Labels") ) ) stop("'d' must have labels.","\n") else sn.d <- attr(d, "Labels")
    if (missing(a) ){
      ab.names <- list("Community1", sn.d)
      a <- matrix(1, 1, n, dimnames = ab.names)
    }
    com <- nrow(a)
    if (!is.matrix(a)) stop("'a' must be a matrix.")
    if (ncol(a) != n )
      stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
    if (is.null(colnames(a) ) ) stop("'a' must have column names", "\n") else sn.a <- colnames(a)
    # check if species labels in 'd' and 'a' match
    if (any(sn.d != sn.a) ) stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).","\n")
    # replace NA in 'a' by 0
    a[which(is.na(a) )] <- 0
    # check if one community has total abundance of zero (no species)
    abun.sum <- apply(a, 1, sum)
    if (any(abun.sum == 0) ) stop("At least one community has zero-sum abundances (no species).","\n")
    # check if one species has total abundance of zero
    abun.sum2 <- apply(a, 2, sum)
    if (any(abun.sum2 == 0) ) stop("At least one species does not occur in any community (zero total abundance across all communities).","\n")
    if (any(is.na(d) ) ) stop("NA's in the distance matrix.","\n")
    A <- matrix(0, ncol = n, nrow = n)
    
    A[row(A) > col(A)] <- -0.5 * d^2
    
    A <- A + t(A)
    
    # gower's double-centering
    G <- bicenter.wt(A)
    e <- eigen(G, symmetric = TRUE)
    
    vectors <- e$vectors
    
    eig <- e$values
    
    # check if 'd' is Euclidean or not
    w0 <- eig[n] / eig[1]
    if (w0 > -tol) r <- sum(eig > (eig[1] * tol)) else r <- length(eig)
    vectors <- vectors[, seq_len(r), drop = FALSE] %*% diag(sqrt(abs(eig <- eig[seq_len(r)])), r)
    
    dimnames(vectors) <- list(colnames(a), NULL)
    pos <- eig > 0
    
    avg.dist.cent <- rep(NA, nrow(a)) ; names(avg.dist.cent) <- row.names(a)
    for (i in seq_len(com)){
      pres <- which(a[i ,] > 0)
      nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F] ) ) )
      if (nb.sp >= 2){
        w <- a[i, pres]
        centroid <- apply(vec, 2, weighted.mean, w = w)
        dist.pos <- sweep(vec[, pos , drop = F], 2, centroid[pos])
        dist.pos <- rowSums(dist.pos^2)
        if (any(!pos) ){
          dist.neg <- sweep(vec[, !pos , drop = F], 2, centroid[!pos])
          dist.neg <- rowSums(dist.neg^2)
        }
        else dist.neg <- 0
        zij <- sqrt(abs(dist.pos - dist.neg) )
        avg.dist.cent[i] <- weighted.mean(zij, w)
      }
      else avg.dist.cent[i] <- 0
    }
    return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
    
  }
  
  
  if (input == "PCoA") {
    
    if ( all( (colnames(a) %in% rownames(d[[1]]) ) != TRUE) )
      stop("the PCoA needs to include all species in the abundance matrix")
    
    a[which(is.na(a) )] <- 0
    # check if one community has total abundance of zero (no species)
    abun.sum <- apply(a, 1, sum)
    if (any(abun.sum == 0) ) stop("At least one community has zero-sum abundances (no species).","\n")
    # check if one species has total abundance of zero
    abun.sum2 <- apply(a, 2, sum)
    if (any(abun.sum2 == 0) ) stop("At least one species does not occur in any community (zero total abundance across all communities).","\n")
    
    com <- nrow(a)
    
    avg.dist.cent <- rep(NA, nrow(a)) ; names(avg.dist.cent) <- row.names(a)
    for (i in seq_len(com)){
      #pres <- which(a[i ,] > 0)
      #nb.sp <- nrow((unique(vec <- d$vectors[pres, , drop = F] ) ) )
      pres <- which(a[i ,] > 0)
      species_com_i <- names(a[i,a[i ,] > 0])
      
      nb.sp <- nrow((unique(vec <- d$vectors[rownames(d$vectors) %in% species_com_i, , drop = F] ) ) )
      
      if (nb.sp >= 2){
        w <- a[i, pres]
        centroid <- apply(vec, 2, weighted.mean, w = w)
        dist.pos <- sweep(vec[, d$pos , drop = F], 2, centroid[d$pos])
        dist.pos <- rowSums(dist.pos^2)
        if (any(!d$pos) ){
          dist.neg <- sweep(vec[, !d$pos , drop = F], 2, centroid[!d$pos])
          dist.neg <- rowSums(dist.neg^2)
        }
        else dist.neg <- 0
        zij <- sqrt(abs(dist.pos - dist.neg) )
        avg.dist.cent[i] <- weighted.mean(zij, w)
      }
      else avg.dist.cent[i] <- 0
    }
    
    return(list(FDis = avg.dist.cent))
    
  }
}


fdisp_updated <- function(d, 
                          a,
                          allow_missing_species_com_input = FALSE,
                          allow_com_func_space_species_discordance = FALSE) {
  
  ### check inputs ###
  if ( !all( (colnames(a) %in% rownames(d[[1]]) ) ) )
    stop("the PCoA needs to include all species in the abundance matrix", call. = FALSE)
  
  if (!allow_com_func_space_species_discordance) {
    if ( !all( (rownames(d[[1]]) %in% colnames(a) ) ) )
      stop("because allow_com_func_space_species_discordance = FALSE, the PCoA cannot include any species that aren't in the community dataset", call. = FALSE)
  }
  
  a[which(is.na(a) )] <- 0
  
  if (any(apply(a, 1, sum) == 0) ) 
    stop("At least one community has zero-sum abundances (no species).","\n", call. = FALSE)
  
  if (!allow_missing_species_com_input && any(apply(a, 2, sum) == 0) ) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).","\n", call. = FALSE)
  
  return(
    apply(a, MARGIN = 1, function(com, d) {
      
      pres <- which(com > 0)
      species_com_i <- names(com[pres])
      
      func_space <- d$vectors[rownames(d$vectors) %in% species_com_i, , drop = FALSE]
      
      #if number of species is at least 2
      if (nrow(func_space) >= 2) {
        w <- com[pres]
        w_reorder <- w[match(rownames(func_space), names(w))]
        
        if (!identical(rownames(func_space), names(w_reorder) )) {
          stop("the order of taxa in the rows of the func space and w_reorder (weights for weighted mean) don't match", call. = FALSE)
        }
        
        centroid <- apply(func_space, 2, weighted.mean, w = w_reorder)
        
        #subtract centroid from each column
        dist.pos <- sweep(func_space[, d$pos, drop = FALSE], 2, STAT = centroid[d$pos], FUN = "-")
        
        #update dist.pos by summing each squared dist.pos
        dist.pos <- rowSums(dist.pos^2)
        
        if (any(!d$pos) ){
          dist.neg <- sweep(func_space[, !d$pos, drop = FALSE], 2, STAT = centroid[!d$pos], FUN = "-")
          dist.neg <- rowSums(dist.neg^2)
        } else {
          dist.neg <- 0
        }
        
        zij <- sqrt(abs(dist.pos - dist.neg) )
        
        if (!identical(names(zij), names(w_reorder) )) {
          stop("the order of taxa in zij and w_reorder (weights for weighted mean) don't match", call. = FALSE)
        }
        
        return(weighted.mean(zij, w_reorder))
      } else {
        return(0)
      }
    }, d = d)
  )
}



calc_fric <- function(funcspace, default_val = NA, options = "Tv", ...) {
  #convex_hull <- NA
  #try(convex_hull <- convhulln(funcspace, options = options, ...), silent = TRUE)
  
  convex_hull <- tryCatch(convhulln(funcspace, options = options, ...),
                          error = function(condition) NA)
  
  if (identical(NA, convex_hull)) {
    return(default_val)
  } else {
    invisible(convex_hull)
  }
}


convhull_vertices <- function(funcspace, default_val = NA) {
  #from the geometry documentation (https://cran.r-project.org/web/packages/geometry/geometry.pdf):
  #if options includes FA, the output will include ...
  #p: The points passed to convnhulln (Alex's commentary: this is just the input functional space)
  #hull: The convex hull, represented as a matrix indexing p, as described above
  
  #Thus, to identify the indices of the points forming the convex hull, we need to identify all the unique
  #values in the hull output. We can then (outside of the convhull_vertices function) identify the actual
  #points by extracting the rows of the functional space (PCoA) using these indices. These results (when
  #sorted in an ascending fashion) perfectly match the vertex identification approach used by FD's dbFD
  #and avoids the need to export and then read in a file.
  
  #convex_hull <- NA
  #try(convex_hull <- convhulln(funcspace, options = 'FA'), silent = TRUE)
  
  convex_hull <- tryCatch(convhulln(funcspace, options = 'FA'),
                          error = function(condition) default_val)
  
  #(1) collapse matrix of indices to a vector
  #(2) extract unique index values
  #(3) sort from smallest --> largest (to match the FD package)
  #if (any(is.na(convex_hull))) {
  if (identical(default_val, convex_hull)) {
    return(default_val)
  } else {
    return(sort(unique(c(convex_hull$hull))) )
  }
}



calc_fdiv <- function(funcspace, com_data, error_val = NA) {
  
  vertices <- convhull_vertices(funcspace, default_val = NA)
  
  if (all(!is.na(vertices))) {
    
    trvertices <- funcspace[vertices, ]
    
    baryv <- apply(trvertices, 2, mean) #center of gravity (mean of vertices of convex full)
    
    #euclidian dstances to Gv (dB)
    
    #for loop based approach to calculating distbaryv (this is the approach used in the FD package)
    #distbaryv <- rep(NA, com_j_richness)
    #for (n in 1:length(com_j)) {
    #distbaryv[n] <- (sum((x.pco_i[n, ] - baryv)^2) )^0.5
    #}
    
    #calculate distbaryv with an apply function (faster than FD's for loop approach)
    distbaryv <- apply(funcspace, MARGIN = 1, FUN = function(x) (sum((x - baryv)^2))^0.5)
    
    #mean of dB values
    meandB <- mean(distbaryv)
    
    #deviations to mean of db
    #devdB <- distbaryv - meandB
    
    #abdev2 and ababsdev2 are weighted by each species' relative abundances
    abundrel <- setNames(com_data[, 2]/sum(com_data[, 2]), com_data[, 1])
    
    #deviations to mean of db
    devdB <- distbaryv - meandB
    devdB <- devdB[match(names(abundrel), names(devdB))]
    
    #message(com_info[positive_eval, 1])
    #message(names(devdB))
    
    if (!identical(names(abundrel), names(devdB))) 
      stop('abundance information does not line up with mean deviation vector')
    
    # relative abundances-weighted mean deviation
    abdev2 <- abundrel*devdB
    
    # relative abundances-weighted mean of absolute deviations
    ababsdev2 <- abundrel*abs(devdB)
    
    #Computation of FDiv 
    return((sum(abdev2) + meandB)/(sum(ababsdev2) + meandB))
    
  } else {
    return(error_val)
  }
}




calc_feve <- function(funcspace, com_data, mst_package = c('vegan', 'ape') ) {
  
  mst_package <- match.arg(mst_package, several.ok = FALSE)
  
  com_j_richness <- nrow(com_data) #species richness of community
  abundrel <- setNames(com_data[, 2]/sum(com_data[, 2]), com_data[, 1])
  
  #calculate distance on the functional space that is ordered to match the order of species in com_data
  tr.dist <- dist(funcspace[match(com_data[, 1], rownames(funcspace)),])
  
  #vegan's spantree (+ custom function) is much faster than ape's mst
  if (mst_package == 'vegan') {
    #vegan package + custom function to calculate/process minimum spanning tree using vegan package's spantree function
    intermediate_mat <- spantree2adjacency_mat(spantree_object = spantree(tr.dist))
    mstvect <- as.dist(intermediate_mat)
  } else {
    #ape's implementation of calculating minimum spanning tree
    linkmst <- mst(tr.dist)
    mstvect <- as.dist(linkmst)
  }
  
  #updated approach to calculating abund2 (faster)
  abund2 <- replicate(com_j_richness, abundrel) + t(replicate(com_j_richness, abundrel))
  dimnames(abund2) <- replicate(2, names(abundrel), simplify = FALSE)
  
  ###FD package's original approach to calculating abund2
  # abund2 <- matrix(0, nrow = com_j_richness, ncol = com_j_richness)
  # dimnames(abund2) <- replicate(2, names(abundrel), simplify = FALSE)
  # for (q in 1:com_j_richness) for (r in 1:com_j_richness) abund2[q, r] <- abundrel[q] + abundrel[r]
  
  abund2vect <- as.dist(abund2)
  
  #check to ensure that all distance matrices have the species in the same order
  if(!all(sapply(list(labels(tr.dist), labels(abund2vect)), FUN = identical, labels(mstvect)))) {
    stop("order of species in the tr.dist, abdund2vec, and mstvect matrices don't match")
  }
  
  EW <- rep(0, com_j_richness - 1)
  flag <- 1
  for (m in seq_len(((com_j_richness - 1)*com_j_richness / 2))) {
    if (mstvect[m] != 0) {
      EW[flag] <- tr.dist[m] / (abund2vect[m])
      flag <- flag + 1
    }
  }
  
  ### computation of the PEW and comparison of PEW with 1/(S - 1) to identify minPEW ###
  #updated version of calculating minPEW
  OdSmO <- 1/(com_j_richness - 1)
  EW_std <- EW/sum(EW)
  minPEW <- sapply(EW_std, function(ew, OdSmO) {
    min(ew, OdSmO)
  }, OdSmO = OdSmO)
  
  #FD package's approach to calculating minPEW
  # minPEW <- rep(0, com_j_richness - 1)
  # OdSmO <- 1 / (com_j_richness - 1)
  # for (l in 1:(com_j_richness - 1) ) {
  #   minPEW[l] <- min((EW[l] / sum(EW)), OdSmO)
  # }
  
  ## premature stop to do some debugging of the FEve code
  #message('community: ', colnames(com_j)[2], '; tr.dist = ', as.character(sum(tr.dist)), '; tr.dist = ', as.character(sum(abund2vect)), '; ', 'EW = ', as.character(sum(EW)), '; FEve: ', ((sum(minPEW)) - OdSmO) / (1 - OdSmO))
  #return(list(EW = EW, mstvect = mstvect, tr.dist = tr.dist, abund2vect = abund2vect))
  
  #computation of FEve
  return(((sum(minPEW)) - OdSmO) / (1 - OdSmO))
}


spantree2adjacency_mat <- function(spantree_object) {
  ### INFORMATION FOR CONSTRUCTING ADJACENCY MATRICES ###
  #(1) examples of adjacency matrices from graphs: https://mathworld.wolfram.com/AdjacencyMatrix.html
  #(2) information on spantree output: (1/4)
  #    "The child node of the parent, starting from parent number two. If there is no link from (2/4) 
  #     the parent, value will be NA and tree is disconnected at the node." (3/4)
  #     source: https://rdrr.io/rforge/vegan/man/spantree.html (4/4)
  
  adjacency_mat <- matrix(data = 0, 
                          nrow = spantree_object$n, ncol = spantree_object$n, 
                          dimnames = list(spantree_object$labels, spantree_object$labels))
  
  for (i in 2:spantree_object$n) {
    adjacency_mat[spantree_object$kid[i - 1], i] <- 1
    adjacency_mat[i, spantree_object$kid[i - 1]] <- 1
  }
  
  return(adjacency_mat)
}


subset_dist <- function(dist_obj, id) {
  dist_mat <- as.matrix(dist_obj)
  return(stats::as.dist(dist_mat[rownames(dist_mat) %in% id, colnames(dist_mat) %in% id]))
}


calc_fd_contr <- function(data, 
                          species = 'all',
                          metrics = c("FRic", "FDis", "FDiv", "FEve"),
                          ntraits,
                          std_traits = TRUE,
                          pcoa_dims,
                          dist_method = c("euclidean", "gowdis", "gawdis"),
                          dist_arg_list = NULL,
                          correction = c("quasieuclid", "cailliez", "lingoes", "sqrt"),
                          feve_calc_type = c("reduced_PCoA", "FD"),
                          mst_package = c('vegan', 'ape'),
                          w_abund = FALSE,
                          calc_contr = TRUE,
                          show_messages = TRUE,
                          include_absent_species = FALSE,
                          unique_id = TRUE,
                          output_PCoA = TRUE) {
  
  #===================================================
  #== Preparatory steps for calculating FD metrics ===
  #===================================================
  
  ### STEP 1. CHECK/PROCESS INPUTS ###
  
  #Checking input arguments
  
  #data
  if (!is.data.frame(data))
    stop("data needs to be a dataframe", call. = FALSE)
  
  rownames(data) <- data[,1] #add species names as the rownames
  
  #species
  if (identical(tolower(species), "all") ) {
    species <- data[,1]
    removed_species <- c()
  } else {
    if (FALSE %in% (species %in% data[,1]) ) {
      stop(paste0("The following species were included in the species argument but are not found in the data: ", paste(species[!(species %in% data[,1])], collapse = ", ") ), call. = FALSE)
    }
    removed_species <- as.character(data[,1][!(data[,1] %in% species) ])
  }
  
  data_preprocessed <- data[(rownames(data) %in% species),] #subset the data down to the specified species
  
  #metrics
  metrics <- match.arg(arg = metrics, several.ok = TRUE)
  
  #ntraits
  if(!is.numeric(ntraits) | (length(ntraits) != 1) )
    stop("ntraits needs to be a single number that indicates the number of traits in the dataframe", call. = FALSE)
  
  #std_traits
  if (!is.logical(std_traits) | length(std_traits) != 1)
    stop("std_traits needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #pcoa_dims
  if(!is.numeric(pcoa_dims) | (length(pcoa_dims) != 1) )
    stop("pcoa_dims needs to be a single number that indicates the number of PCoA axes to use in the functional space", call. = FALSE)
  
  #dist_method
  dist_method <- match.arg(arg = dist_method, several.ok = FALSE)
  
  #https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame
  if ((dist_method == "euclidean") & ( !all(sapply(data[, 2:(ntraits + 1)], is.numeric)) | !all(!is.na(data[, 2:(ntraits + 1)])) ) ) {
    dist_method <- "gowdis"
    
    if (show_messages) {
      message("Since not all traits are numeric, distance matrix calculation will use Gower's distance instead of Euclidean distance (the method that was specified)")
    }
  }
  
  #correction
  correction <- match.arg(arg = correction, several.ok = FALSE)
  
  #feve_calc_type
  feve_calc_type <- match.arg(arg = feve_calc_type, several.ok = FALSE)
  
  #mst_package
  mst_package <- match.arg(mst_package, several.ok = FALSE)
  
  #w_abund
  if (!is.logical(w_abund) | length(w_abund) != 1 )
    stop("w_abund needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  #calc_contr
  if (!is.logical(calc_contr) | length(calc_contr) != 1 )
    stop("calc_contr needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #show_messages
  if (!is.logical(show_messages) | length(show_messages) != 1 )
    stop("show_messages needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #include_absent_species
  if (!is.logical(include_absent_species) | length(include_absent_species) != 1 )
    stop("include_absent_species needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  #unique_id
  if (!is.logical(unique_id) | length(unique_id) != 1 )
    stop("unique_id needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE) 
  
  #output_PCoA
  if (!is.logical(output_PCoA) | length(output_PCoA) != 1 )
    stop("output_PCoA needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  
  ### Identify and handle communities with insufficient species richness ###
  species_rich_dataframe <- as.data.frame(apply(data_preprocessed[, -c(1:(1 + ntraits)), drop = FALSE], MARGIN = 2, function(x) { length(x[x > 0 ]) } ))
  species_rich_dataframe$community <- rownames(species_rich_dataframe)
  species_rich_dataframe <- species_rich_dataframe[, c(2, 1)]
  colnames(species_rich_dataframe) <- c("community", "species_richness")
  
  if (isTRUE(calc_contr)) {
    #the FD contribution metrics require at least TWO more PCoA dimensions than the minimum community species richness
    insufficient_coms <- species_rich_dataframe[which(species_rich_dataframe$species_richness < (pcoa_dims + 2) ),]$community
  } else {
    #the community-level FD metrics (without calculating contribution) require at least ONE more PCoA dimension than the minimum community species richness
    insufficient_coms <- species_rich_dataframe[which(species_rich_dataframe$species_richness < (pcoa_dims + 1) ),]$community
  }
  
  ifelse(length(insufficient_coms) == 0,
         insufficient_com_info <- "All communities had sufficient species richness based on the number of specified PCoA axes",
         insufficient_com_info <- species_rich_dataframe[species_rich_dataframe$community %in% insufficient_coms,])
  
  #remove communities with insufficient species richness
  data_processed <- data_preprocessed[,!(colnames(data_preprocessed) %in% insufficient_coms)]
  
  
  ### Identify and handle species that are not found in any communities ###
  
  #find species that are absent across all communities
  absent_species <- as.character(data_processed[apply( data_processed[,-c(1:(1 + ntraits)), drop = FALSE], MARGIN = 1, function(x) all(x == 0) ), 1])
  
  if (include_absent_species) {
    #if include_absent_species == TRUE, report the absent species but don't remove them from the dataset 
    ifelse( (length(absent_species) == 0) & (length(removed_species) == 0),
            removed_species_info <- paste0("No species were absent from all communities or were manually excluded"),
            removed_species_info <- list(message = "The absent species were not removed because include_absent_species == TRUE",
                                         removed_species = removed_species,
                                         missing_species = absent_species))
  } else {
    #if include_absent_species == FALSE, report the absent species and remove them from the dataset
    
    #remove absent species 
    data_processed <- data_processed[ !(data_processed[, 1] %in% absent_species) ,]
    
    #absent species message
    ifelse( (length(absent_species) == 0) & (length(removed_species) == 0) ,
            removed_species_info <- paste0("No species were absent from all communities or were manually excluded"),
            removed_species_info <- list(message = "The absent species were removed because include_absent_species == FALSE",
                                         removed_species = removed_species,
                                         absent_species = absent_species) )
  }
  
  if (isTRUE(show_messages)) {
    ifelse(length(data[,1][!(data[, 1] %in% species) ]) != 0,
           removed_spp_print <- paste(data[, 1][!(data[, 1] %in% species) ], collapse = ", "),
           removed_spp_print <- "no species were manually excluded")
    
    ifelse(length(absent_species) != 0,
           absent_species_print <- paste(absent_species, collapse = ", "),
           absent_species_print <- "no species were absent across all communities")
    
    #including a line break in a printed message
    #https://stackoverflow.com/questions/18247990/print-multiple-lines-without-printing-the-print
    message(paste0("### REMOVED THE FOLLOWING SPECIES: ### \nManually specified: ", removed_spp_print,
                   "; \nAbsent across all communities: ", absent_species_print, "\n" ) )
  }
  
  
  ### Step 2: CREATE FUNCTIONAL SPACE ###
  
  #2a: Process trait data
  if (isTRUE(std_traits) & dist_method == "euclidean") {
    #if std_traits and distance matrix construction is euclidean
    trait_data <- apply(data_processed[, 2:(ntraits + 1), drop = FALSE], MARGIN = 2, function(x) standardize_func(x))
    
  } else {
    
    #if std_traits != TRUE and/or distance matrix construction is gower
    trait_data <- data_processed[, 2:(ntraits + 1), drop = FALSE]
  }
  
  
  #2b: Calculate trait dissimilarity matrix
  dist_mat <- dist_wrapper(df = trait_data,
                           method = switch(dist_method, euclidean = 'dist', gowdis = 'gowdis', gawdis = 'gawdis'),
                           arg_list = dist_arg_list)
  
  #distance_matrix_list <- list(dist_mat = dist_mat$dist,
  #                             dist_mat_info = dist_mat$function_info)
  
  distance_matrix_list <- list(dist_mat = dist_mat,
                               dist_mat_info = dist_method)
  
  
  #2c: Determine if correction is needed, apply correction, compute PCoA
  
  if (is.euclid(distance_matrix_list[["dist_mat"]])) {
    correction_message <- "No correction applied"
    euclid_distmat <- distance_matrix_list[["dist_mat"]]
    
  } else {
    correction_message <- paste0("Correction applied: ", correction)
    
    euclid_distmat <- switch(correction,
                             quasieuclid = {quasieuclid(distance_matrix_list[["dist_mat"]])},
                             cailliez = {cailliez(distance_matrix_list[["dist_mat"]])},
                             lingoes = {lingoes(distance_matrix_list[["dist_mat"]])},
                             sqrt = {sqrt(distance_matrix_list[["dist_mat"]])})
  }
  
  x.pco_init <- ade4::dudi.pco(euclid_distmat, scannf = FALSE, full = TRUE)
  
  
  if (pcoa_dims > ncol(x.pco_init$li)) {
    #THIS NEEDS TO BE FIXED. CHANGE TO THE MAX NUMBER OF PCOA AXES THAT IS ACCEPTABLE BASED ON SPECIES RICHNESS OF THE ASSEMBLAGE
    x.pco_full <- x.pco_init$li
    
    if (isTRUE(show_messages)) {
      message("The number of specified PCoA axes is greater than the number of available PCoA axes. Using the all available PCoA axes in the functional space.")
    }
    
  } else {
    x.pco_full <- x.pco_init$li[, seq_len(pcoa_dims), drop = FALSE]
  }
  
  #if the FD calculations are only going to consider presence/absence information (w_abund == FALSE), change all non-zero abundances to 1
  if (isFALSE(w_abund)) {
    data_processed[, (ntraits + 2):ncol(data_processed)][data_processed[,(ntraits + 2):ncol(data_processed)] != 0] <- 1
  }
  
  
  #==================================================
  #== Calculation of functional diversity metrics ===
  #==================================================
  
  if ("FDis" %in% metrics) {
    # Creating input files for FDis calculations
    #FOLLOWING THE FD PACKAGE, I AM USING THE UNCORRECTED DISTANCE MATRIX SINCE THE FDIS PREPARATION INCLUDES A CORRECTION
    #fdis_input <- FDis_PCoA_input(d = euclid_distmat)
    fdis_input <- FDis_PCoA_input(d = distance_matrix_list[["dist_mat"]])
  }
  
  com_diversity_list <- list() #list to store the FD metrics for each community
  
  if (isTRUE(show_messages)) {
    message("starting functional diversity metric calculations")
    progress_bar <- txtProgressBar(min = 0, max = (ncol(data_processed) - (1 + ntraits) ), style = 3, char = "*") #initiate progress bar
  }
  
  #calculate functional diversity metrics
  for (j in seq_len( ncol(data_processed) - (ntraits + 1)) ) {
    
    starting_val <- 1 + ntraits #starting position of community data in the dataframe
    
    #subset dataframe down species found in community j and the species and jth community columns
    com_j <- data_processed[which(data_processed[, (starting_val + j)] > 0), c(1, (starting_val + j) )]
    
    species_com_j <- com_j[, 1] #species found in community j
    rownames(com_j) <- species_com_j #add species to the rownames
    com_j_richness <- nrow(com_j) #species richness of community
    
    if (isTRUE(calc_contr)) {
      diversity_data <- data.frame(species = species_com_j)
      diversity_data$community <- colnames(data_processed)[starting_val + j]
      diversity_data$species_richness <- com_j_richness
      
    } else {
      diversity_data <- data.frame(community = colnames(data_processed)[starting_val + j])
      diversity_data$species_richness <- com_j_richness
    }
    
    
    #subset the full PCoA matrix down to the species in community j
    x.pco_iterativedrop_initial <- x.pco_full[which(rownames(x.pco_full) %in% species_com_j),, drop = FALSE]
    x.pco_iterativedrop <- x.pco_iterativedrop_initial[match(species_com_j, rownames(x.pco_iterativedrop_initial)),, drop = FALSE]
    
    if ("FRic" %in% metrics) {
      ### full FRic  (includes all species in community) ###
      
      #diversity_data$full_FRic <- convhulln(x.pco_iterativedrop, options = c("QJ", "FA"))$vol
      full_fric <- calc_fric(funcspace = x.pco_iterativedrop, default_val = NA, options = c("QJ", "FA"))
      diversity_data$full_FRic <- ifelse(any(is.na(full_fric)), 9999, full_fric$vol)
      
      if (isTRUE(calc_contr)) {
        ### each community constituent's contribution to FRic ###
        for (i in seq_len(length(species_com_j))) {
          exclude <- species_com_j[i] #focal species to exclude from FRic calculation
          x.pco_subset_alt <- x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE]
          
          #FRic_reduced <- convhulln(x.pco_subset_alt, options = c("QJ", "FA"))$vol
          #diversity_data[i,'partial_FRic'] <- FRic_reduced
          
          partial_fric <- calc_fric(funcspace = x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE], 
                                    default_val = NA, 
                                    options = c("QJ", "FA"))
          diversity_data[i,'partial_FRic'] <- ifelse(any(is.na(partial_fric)), 9999, partial_fric$vol)
        }
      }
    }
    
    
    if ("FDis" %in% metrics) {
      ### Full FDis  (includes all species in community) ###
      
      matrix_FDis <- t(as.matrix(com_j[,2]))
      
      colnames(matrix_FDis) <- species_com_j
      
      #comj_species changed to as.character(com_j$species) because when only one community is simulated, the com_j$species vector is a factor, which messes up the dis_subset function (it subsets the levels not the actual values). This inexplicably is not a problem when more than one community is simulated (it is a character vector)
      #diversity_data$full_FDis_V1 <- fdisp(dist_subset(euclid_distmat, as.character(species_com_j) ), matrix_FDis)$FDis
      
      #diversity_data$full_FDis <- fdisp_alt(d = fdis_input, a = matrix_FDis, input = "PCoA", tol = 1e-07)$FDis
      diversity_data$full_FDis <- fdisp_updated(d = fdis_input, 
                                                a = matrix_FDis,
                                                allow_missing_species_com_input = TRUE,
                                                allow_com_func_space_species_discordance = TRUE)
      
      if (isTRUE(calc_contr)) {
        ### each community constituent's contribution to FDis ###
        for (i in seq_len(length(species_com_j))) {
          exclude <- species_com_j[i] #focal species to exclude from FDis calculation
          com_j_names_reduced <- species_com_j[!(species_com_j %in% exclude)]
          
          matrix_FDis_reduced <- t(as.matrix(com_j[ !(com_j[,1] %in% exclude), 2]))
          colnames(matrix_FDis_reduced) <- com_j_names_reduced
          
          #diversity_data[i,'partial_FDis'] <- fdisp_alt(d = fdis_input, a = matrix_FDis_reduced, input = "PCoA", tol = 1e-07)$FDis
          diversity_data[i,'partial_FDis'] <- fdisp_updated(d = fdis_input, 
                                                            a = matrix_FDis_reduced,
                                                            allow_missing_species_com_input = TRUE,
                                                            allow_com_func_space_species_discordance = TRUE)
        }
      }
    }
    
    
    if ("FDiv" %in% metrics) {
      ### Full FDiv (includes all species in community) ###
      
      diversity_data$full_FDiv <- calc_fdiv(funcspace = x.pco_iterativedrop, 
                                            com_data = com_j, 
                                            error_val = 9999)
      
      if (isTRUE(calc_contr)) {
        
        ### Each community constituent's contribution to FDiv ###
        for (i in seq_len(length(species_com_j))) {
          ### FDiv reduction (reduced_FDiv) --> fully recalculates FDiv after dropping each species ###
          #species contributions calculated on reduced_FDiv represent each species' contribution to the observed FDiv of the community
          
          exclude <- species_com_j[i] #focal species to exclude from FDiv calculation
          x.pco_subset_alt <- x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE]
          
          diversity_data[i,'partial_FDiv'] <- calc_fdiv(funcspace = x.pco_subset_alt, 
                                                        com_data = com_j[!(com_j[, 1] %in% exclude),, drop = FALSE], 
                                                        error_val = 9999)
        }
      }
    }
    
    
    if ("FEve" %in% metrics) {
      ### Full FEve (includes all species in community) ###
      
      if (feve_calc_type == "reduced_PCoA") {
        #alternative approach to computing the minimum spanning tree that uses the same number of PCoA axes as the FRic and FDiv calculations
        diversity_data$full_FEve <- calc_feve(funcspace = x.pco_iterativedrop, com_data = com_j, mst_package = mst_package)
        
      } else {
        ### THE FD PACKAGE USES ALL AXES OF THE PCoA TO COMPUTE THE MINIMUM SPANNING TREE ###
        #to replicate the FD package method, uncomment these two lines and comment out the alternative tr.dist computation directly below
        FEve_PCoA <- x.pco_init$li[which(rownames(x.pco_init$li) %in% species_com_j),, drop = FALSE]
        #FEve_PCoA <- FEve_PCoA[match(species_com_j, rownames(FEve_PCoA)),]
        
        diversity_data$full_FEve <- calc_feve(funcspace = FEve_PCoA, com_data = com_j, mst_package = mst_package)
      }
      
      if (isTRUE(calc_contr)) {
        
        for (i in seq_len(length(species_com_j))) {
          
          ### each community constituent's contribution to FEve ###
          exclude <- species_com_j[i] #focal species to exclude from FEve calculation
          
          if (feve_calc_type == "reduced_PCoA") {
            #alternative approach to computing the minimum spanning tree that uses the same PCoA axes as the FRic and FDiv calculations
            diversity_data[i,'partial_FEve'] <- calc_feve(funcspace = x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE], 
                                                          com_data = com_j[ !(species_com_j %in% exclude),, drop = FALSE], 
                                                          mst_package = mst_package)
            
          } else {
            #THE APPROACH USED IN THE FD PACKAGE:
            ### THE FD PACKAGE USES ALL AXES OF THE PCoA TO COMPUTE THE MINIMUM SPANNING TREE ###
            diversity_data[i,'partial_FEve'] <- calc_feve(funcspace = FEve_PCoA[!(row.names(FEve_PCoA) %in% exclude),, drop = FALSE], 
                                                          com_data = com_j[ !(species_com_j %in% exclude),, drop = FALSE], 
                                                          mst_package = mst_package)
          }
        }
      }
    }
    
    
    ### Calculate the unstandardized and standardized functional contribution metrics ###
    if (isTRUE(calc_contr)) {
      for (m in metrics) {
        diversity_data[, paste0(m, "_prop")] <-  diversity_data[, paste0("partial_", m)]/diversity_data[, paste0("full_", m)]
        diversity_data[, paste0(m, "_contr")] <-  1 - diversity_data[, paste0(m, "_prop")]
        diversity_data[, paste0(m, "_contr_std")] <- standardize_func(diversity_data[, paste0(m, "_contr")])
      }
    }
    
    com_diversity_list[[j]] <- diversity_data #add the data/calculations for community j to com_diversity_list
    
    if (isTRUE(show_messages)) setTxtProgressBar(progress_bar, j) #update progress bar
    
  }
  
  
  #==========================================
  #== Preparation of materials for output ===
  #==========================================
  
  #collapse the dataframes within com_diversity_list into a single dataframe
  com_diversity_list_collapsed <- do.call(rbind, com_diversity_list)
  
  #add unique identifier (unique_id) if unique_id == TRUE
  if (isTRUE(calc_contr) & isTRUE(unique_id)) com_diversity_list_collapsed$unique_id <- paste(com_diversity_list_collapsed$species, com_diversity_list_collapsed$community, sep = "_")
  
  #identify the problem communities whose convex hull couldn't be formed (full_FDiv, full_FRic, partial_FDiv, and/or partial_FRic = 9999)
  included_cols <- c('full_FDiv', 'partial_FDiv', 'full_FRic', 'partial_FRic')[c('full_FDiv', 'partial_FDiv', 'full_FRic', 'partial_FRic') %in% colnames(com_diversity_list_collapsed)]
  if (length(included_cols) > 0) {
    problem_coms <- com_diversity_list_collapsed[apply(com_diversity_list_collapsed[, included_cols, drop = FALSE], 1, function(x) 9999 %in% x),]$community
  } else {
    problem_coms <- vector()
  }
  
  #remove problem communities
  com_diversity_list_FINAL <- com_diversity_list_collapsed[!(com_diversity_list_collapsed$community %in% problem_coms),, drop = FALSE]
  
  #report problem communities
  ifelse(length(problem_coms) == 0, 
         problem_com_info <- "No problem communities",
         problem_com_info <- problem_coms)
  
  #export functional diversity calculations dataframe and auxiliary information as a list
  if (isTRUE(output_PCoA)) {
    #include PCoA in output
    return(list(FD_dataframe = com_diversity_list_FINAL,
                distance_matrix_info = distance_matrix_list,
                correction = correction_message,
                communities_with_insufficient_richness = insufficient_com_info,
                absent_species_info = removed_species_info,
                failed_communities = problem_com_info,
                PCoA = x.pco_full))    
    
  } else {
    #don't include PCoA in output
    return(list(FD_dataframe = com_diversity_list_FINAL,
                distance_matrix_info = distance_matrix_list,
                correction = correction_message,
                communities_with_insufficient_richness = insufficient_com_info,
                absent_species_info = removed_species_info,
                failed_communities = problem_com_info))
  }
  
} #end of calc_fd_contr function




calc_fd_contr_inputdist <- function(data,
                                    trait_dist,
                                    species = 'all',
                                    metrics = c("FRic", "FDis", "FDiv", "FEve"),
                                    std_traits = TRUE,
                                    pcoa_dims,
                                    dist_arg_list = NULL,
                                    correction = c("quasieuclid", "cailliez", "lingoes", "sqrt"),
                                    feve_calc_type = c("reduced_PCoA", "FD"),
                                    mst_package = c('vegan', 'ape'),
                                    w_abund = FALSE,
                                    calc_contr = TRUE,
                                    show_messages = TRUE,
                                    include_absent_species = FALSE,
                                    unique_id = TRUE,
                                    output_PCoA = TRUE) {
  
  #===================================================
  #== Preparatory steps for calculating FD metrics ===
  #===================================================
  
  ### STEP 1. CHECK/PROCESS INPUTS ###
  
  #Checking input arguments
  
  #data
  #if (!is.list(data) | all(!sapply(data, is.data.frame)) | all(!sapply(data, function(x) class(x) == 'dist')))
  #  stop("data needs to be a list that contains a dataframe (of community info) and a distance object of trait distances", call. = FALSE)
  if (!is.data.frame(data))
    stop("data needs to be a dataframe", call. = FALSE)
  
  if (class(trait_dist) != 'dist')
    stop("trait_dist needs to be an object of class dist", call. = FALSE)
  
  community_df <- data
  rownames(community_df) <- community_df[, 1] #add species names as the rownames
  
  #trait_dist <- data[[which(sapply(data, function(x) class(x) == 'dist'))]]
  
  #make sure that all taxa in community dataset are found in the distance matrix
  if ( any(!(community_df[, 1] %in% colnames(as.matrix(trait_dist)))))
    stop('there are species in the community dataset that are not in the dissimilarity matrix')
  
  
  #species
  if (identical(tolower(species), "all") ) {
    species <- community_df[,1]
    removed_species <- c()
  } else {
    if (FALSE %in% (species %in% community_df[,1]) ) {
      stop(paste0("The following species were included in the species argument but are not found in the data: ", paste(species[!(species %in% data[,1])], collapse = ", ") ), call. = FALSE)
    }
    removed_species <- as.character(community_df[,1][!(community_df[,1] %in% species) ])
  }
  
  data_preprocessed <- community_df[(rownames(community_df) %in% species),] #subset the data down to the specified species
  dist_subset <- subset_dist(dist_obj = trait_dist, id = data_preprocessed[,1])
  
  ### double check harmonization of species between community dataset and distance matrix ###
  if (!identical(sort(data_preprocessed[,1]), sort(colnames(as.matrix(dist_subset)))))
    stop("The distance matrix and community dataset do contain identical sets of taxa")
  
  #metrics
  metrics <- match.arg(arg = metrics, several.ok = TRUE)
  
  #ntraits
  #if(!is.numeric(ntraits) | (length(ntraits) != 1) )
  #  stop("ntraits needs to be a single number that indicates the number of traits in the dataframe", call. = FALSE)
  
  #std_traits
  #if (!is.logical(std_traits) | length(std_traits) != 1)
  #  stop("std_traits needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #pcoa_dims
  if(!is.numeric(pcoa_dims) | (length(pcoa_dims) != 1) )
    stop("pcoa_dims needs to be a single number that indicates the number of PCoA axes to use in the functional space", call. = FALSE)
  
  #dist_method
  #dist_method <- match.arg(arg = dist_method, several.ok = FALSE)
  
  #https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame
  #if ((dist_method == "euclidean") & ( !all(sapply(data[, 2:(ntraits + 1)], is.numeric)) | !all(!is.na(data[, 2:(ntraits + 1)])) ) ) {
  #  dist_method <- "gowdis"
  #  
  #  if (show_messages) {
  #    message("Since not all traits are numeric, distance matrix calculation will use Gower's distance instead of Euclidean distance (the method that was specified)")
  #  }
  #}
  
  #correction
  correction <- match.arg(arg = correction, several.ok = FALSE)
  
  #feve_calc_type
  feve_calc_type <- match.arg(arg = feve_calc_type, several.ok = FALSE)
  
  #mst_package
  mst_package <- match.arg(mst_package, several.ok = FALSE)
  
  #w_abund
  if (!is.logical(w_abund) | length(w_abund) != 1 )
    stop("w_abund needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  #calc_contr
  if (!is.logical(calc_contr) | length(calc_contr) != 1 )
    stop("calc_contr needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #show_messages
  if (!is.logical(show_messages) | length(show_messages) != 1 )
    stop("show_messages needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE)
  
  #include_absent_species
  if (!is.logical(include_absent_species) | length(include_absent_species) != 1 )
    stop("include_absent_species needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  #unique_id
  if (!is.logical(unique_id) | length(unique_id) != 1 )
    stop("unique_id needs to be a single logical value (TRUE or FALSE). Default is TRUE.", call. = FALSE) 
  
  #output_PCoA
  if (!is.logical(output_PCoA) | length(output_PCoA) != 1 )
    stop("output_PCoA needs to be a single logical value (TRUE or FALSE). Default is FALSE.", call. = FALSE)
  
  
  ### Identify and handle communities with insufficient species richness ###
  species_rich_dataframe <- as.data.frame(apply(data_preprocessed[, -1, drop = FALSE], MARGIN = 2, function(x) { length(x[x > 0 ]) } ))
  species_rich_dataframe$community <- rownames(species_rich_dataframe)
  species_rich_dataframe <- species_rich_dataframe[, c(2, 1)]
  colnames(species_rich_dataframe) <- c("community", "species_richness")
  
  if (isTRUE(calc_contr)) {
    #the FD contribution metrics require at least TWO more PCoA dimensions than the minimum community species richness
    insufficient_coms <- species_rich_dataframe[which(species_rich_dataframe$species_richness < (pcoa_dims + 2) ),]$community
  } else {
    #the community-level FD metrics (without calculating contribution) require at least ONE more PCoA dimension than the minimum community species richness
    insufficient_coms <- species_rich_dataframe[which(species_rich_dataframe$species_richness < (pcoa_dims + 1) ),]$community
  }
  
  ifelse(length(insufficient_coms) == 0,
         insufficient_com_info <- "All communities had sufficient species richness based on the number of specified PCoA axes",
         insufficient_com_info <- species_rich_dataframe[species_rich_dataframe$community %in% insufficient_coms,])
  
  #remove communities with insufficient species richness
  data_processed <- data_preprocessed[,!(colnames(data_preprocessed) %in% insufficient_coms)]
  
  
  ### Identify and handle species that are not found in any communities ###
  
  #find species that are absent across all communities
  absent_species <- as.character(data_processed[apply( data_processed[, -1, drop = FALSE], MARGIN = 1, function(x) all(x == 0) ), 1])
  
  if (isTRUE(include_absent_species)) {
    #if include_absent_species == TRUE, report the absent species but don't remove them from the dataset 
    ifelse( (length(absent_species) == 0) & (length(removed_species) == 0),
            removed_species_info <- paste0("No species were absent from all communities or were manually excluded"),
            removed_species_info <- list(message = "The absent species were not removed because include_absent_species == TRUE",
                                         removed_species = removed_species,
                                         missing_species = absent_species))
  } else {
    #if include_absent_species == FALSE, report the absent species and remove them from the dataset
    
    #remove absent species 
    data_processed <- data_processed[ !(data_processed[, 1] %in% absent_species) ,]
    
    #absent species message
    ifelse( (length(absent_species) == 0) & (length(removed_species) == 0) ,
            removed_species_info <- paste0("No species were absent from all communities or were manually excluded"),
            removed_species_info <- list(message = "The absent species were removed because include_absent_species == FALSE",
                                         removed_species = removed_species,
                                         absent_species = absent_species) )
  }
  
  if (isTRUE(show_messages)) {
    ifelse(length(community_df[,1][!(community_df[, 1] %in% species) ]) != 0,
           removed_spp_print <- paste(community_df[, 1][!(community_df[, 1] %in% species) ], collapse = ", "),
           removed_spp_print <- "no species were manually excluded")
    
    ifelse(length(absent_species) != 0,
           absent_species_print <- paste(absent_species, collapse = ", "),
           absent_species_print <- "no species were absent across all communities")
    
    #including a line break in a printed message
    #https://stackoverflow.com/questions/18247990/print-multiple-lines-without-printing-the-print
    message(paste0("### REMOVED THE FOLLOWING SPECIES: ### \nManually specified: ", removed_spp_print,
                   "; \nAbsent across all communities: ", absent_species_print, "\n" ) )
  }
  
  
  ### Step 2: CREATE FUNCTIONAL SPACE ###
  
  #2a: Process trait data
  # if (isTRUE(std_traits) & dist_method == "euclidean") {
  #   #if std_traits and distance matrix construction is euclidean
  #   trait_data <- apply(data_processed[, 2:(ntraits + 1), drop = FALSE], MARGIN = 2, function(x) standardize_func(x))
  #   
  # } else {
  #   
  #   #if std_traits != TRUE and/or distance matrix construction is gower
  #   trait_data <- data_processed[, 2:(ntraits + 1), drop = FALSE]
  # }
  # 
  
  #2b: Calculate trait dissimilarity matrix
  #dist_mat <- dist_wrapper(df = trait_data,
  #                         method = switch(dist_method, euclidean = 'dist', gowdis = 'gowdis', gawdis = 'gawdis'),
  #                         arg_list = dist_arg_list)
  
  distance_matrix_list <- list(dist_mat = dist_subset,
                               dist_mat_info = 'user input dissimilarity matrix')
  
  
  #2c: Determine if correction is needed, apply correction, compute PCoA
  
  if (is.euclid(distance_matrix_list[["dist_mat"]])) {
    correction_message <- "No correction applied"
    euclid_distmat <- distance_matrix_list[["dist_mat"]]
    
  } else {
    correction_message <- paste0("Correction applied: ", correction)
    
    euclid_distmat <- switch(correction,
                             quasieuclid = {quasieuclid(distance_matrix_list[["dist_mat"]])},
                             cailliez = {cailliez(distance_matrix_list[["dist_mat"]])},
                             lingoes = {lingoes(distance_matrix_list[["dist_mat"]])},
                             sqrt = {sqrt(distance_matrix_list[["dist_mat"]])})
  }
  
  x.pco_init <- ade4::dudi.pco(euclid_distmat, scannf = FALSE, full = TRUE)
  
  
  if (pcoa_dims > ncol(x.pco_init$li)) {
    #THIS NEEDS TO BE FIXED. CHANGE TO THE MAX NUMBER OF PCOA AXES THAT IS ACCEPTABLE BASED ON SPECIES RICHNESS OF THE ASSEMBLAGE
    x.pco_full <- x.pco_init$li
    
    if (isTRUE(show_messages)) {
      message("The number of specified PCoA axes is greater than the number of available PCoA axes. Using the all available PCoA axes in the functional space.")
    }
    
  } else {
    x.pco_full <- x.pco_init$li[, seq_len(pcoa_dims), drop = FALSE]
  }
  
  #if the FD calculations are only going to consider presence/absence information (w_abund == FALSE), change all non-zero abundances to 1
  if (isFALSE(w_abund)) {
    #data_processed[, (ntraits + 2):ncol(data_processed)][data_processed[,(ntraits + 2):ncol(data_processed)] != 0] <- 1
    data_processed[, 2:ncol(data_processed)][data_processed[, 2:ncol(data_processed)] != 0] <- 1
  }
  
  
  #==================================================
  #== Calculation of functional diversity metrics ===
  #==================================================
  
  if ("FDis" %in% metrics) {
    # Creating input files for FDis calculations
    #FOLLOWING THE FD PACKAGE, I AM USING THE UNCORRECTED DISTANCE MATRIX SINCE THE FDIS PREPARATION INCLUDES A CORRECTION
    #fdis_input <- FDis_PCoA_input(d = euclid_distmat)
    fdis_input <- FDis_PCoA_input(d = distance_matrix_list[["dist_mat"]])
  }
  
  com_diversity_list <- list() #list to store the FD metrics for each community
  
  if (isTRUE(show_messages)) {
    message("starting functional diversity metric calculations")
    progress_bar <- txtProgressBar(min = 0, max = (ncol(data_processed) - 1 ), style = 3, char = "*") #initiate progress bar
  }
  
  #calculate functional diversity metrics
  for (j in seq_len( ncol(data_processed) - 1) ) {
    
    starting_val <- 1 #starting position of community data in the dataframe
    
    #subset dataframe down species found in community j and the species and jth community columns
    com_j <- data_processed[which(data_processed[, (starting_val + j)] > 0), c(1, (starting_val + j) )]
    
    species_com_j <- com_j[, 1] #species found in community j
    rownames(com_j) <- species_com_j #add species to the rownames
    com_j_richness <- nrow(com_j) #species richness of community
    
    if (isTRUE(calc_contr)) {
      diversity_data <- data.frame(species = species_com_j)
      diversity_data$community <- colnames(data_processed)[starting_val + j]
      diversity_data$species_richness <- com_j_richness
      
    } else {
      diversity_data <- data.frame(community = colnames(data_processed)[starting_val + j])
      diversity_data$species_richness <- com_j_richness
    }
    
    
    #subset the full PCoA matrix down to the species in community j
    x.pco_iterativedrop_initial <- x.pco_full[which(rownames(x.pco_full) %in% species_com_j),, drop = FALSE]
    x.pco_iterativedrop <- x.pco_iterativedrop_initial[match(species_com_j, rownames(x.pco_iterativedrop_initial)),, drop = FALSE]
    
    if ("FRic" %in% metrics) {
      ### full FRic  (includes all species in community) ###
      
      #diversity_data$full_FRic <- convhulln(x.pco_iterativedrop, options = c("QJ", "FA"))$vol
      full_fric <- calc_fric(funcspace = x.pco_iterativedrop, default_val = NA, options = c("QJ", "FA"))
      diversity_data$full_FRic <- ifelse(any(is.na(full_fric)), 9999, full_fric$vol)
      
      if (isTRUE(calc_contr)) {
        ### each community constituent's contribution to FRic ###
        for (i in seq_len(length(species_com_j))) {
          exclude <- species_com_j[i] #focal species to exclude from FRic calculation
          x.pco_subset_alt <- x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE]
          
          #FRic_reduced <- convhulln(x.pco_subset_alt, options = c("QJ", "FA"))$vol
          #diversity_data[i,'partial_FRic'] <- FRic_reduced
          
          partial_fric <- calc_fric(funcspace = x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE], 
                                    default_val = NA, 
                                    options = c("QJ", "FA"))
          diversity_data[i,'partial_FRic'] <- ifelse(any(is.na(partial_fric)), 9999, partial_fric$vol)
        }
      }
    }
    
    
    if ("FDis" %in% metrics) {
      ### Full FDis  (includes all species in community) ###
      
      matrix_FDis <- t(as.matrix(com_j[,2]))
      
      colnames(matrix_FDis) <- species_com_j
      
      #comj_species changed to as.character(com_j$species) because when only one community is simulated, the com_j$species vector is a factor, which messes up the dis_subset function (it subsets the levels not the actual values). This inexplicably is not a problem when more than one community is simulated (it is a character vector)
      #diversity_data$full_FDis_V1 <- fdisp(dist_subset(euclid_distmat, as.character(species_com_j) ), matrix_FDis)$FDis
      
      #diversity_data$full_FDis <- fdisp_alt(d = fdis_input, a = matrix_FDis, input = "PCoA", tol = 1e-07)$FDis
      diversity_data$full_FDis <- fdisp_updated(d = fdis_input, 
                                                a = matrix_FDis,
                                                allow_missing_species_com_input = TRUE,
                                                allow_com_func_space_species_discordance = TRUE)
      
      if (isTRUE(calc_contr)) {
        ### each community constituent's contribution to FDis ###
        for (i in seq_len(length(species_com_j))) {
          exclude <- species_com_j[i] #focal species to exclude from FDis calculation
          com_j_names_reduced <- species_com_j[!(species_com_j %in% exclude)]
          
          matrix_FDis_reduced <- t(as.matrix(com_j[ !(com_j[,1] %in% exclude), 2]))
          colnames(matrix_FDis_reduced) <- com_j_names_reduced
          
          #diversity_data[i,'partial_FDis'] <- fdisp_alt(d = fdis_input, a = matrix_FDis_reduced, input = "PCoA", tol = 1e-07)$FDis
          diversity_data[i,'partial_FDis'] <- fdisp_updated(d = fdis_input, 
                                                            a = matrix_FDis_reduced,
                                                            allow_missing_species_com_input = TRUE,
                                                            allow_com_func_space_species_discordance = TRUE)
        }
      }
    }
    
    
    if ("FDiv" %in% metrics) {
      ### Full FDiv (includes all species in community) ###
      
      diversity_data$full_FDiv <- calc_fdiv(funcspace = x.pco_iterativedrop, 
                                            com_data = com_j, 
                                            error_val = 9999)
      
      if (isTRUE(calc_contr)) {
        
        ### Each community constituent's contribution to FDiv ###
        for (i in seq_len(length(species_com_j))) {
          ### FDiv reduction (reduced_FDiv) --> fully recalculates FDiv after dropping each species ###
          #species contributions calculated on reduced_FDiv represent each species' contribution to the observed FDiv of the community
          
          exclude <- species_com_j[i] #focal species to exclude from FDiv calculation
          x.pco_subset_alt <- x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE]
          
          diversity_data[i,'partial_FDiv'] <- calc_fdiv(funcspace = x.pco_subset_alt, 
                                                        com_data = com_j[!(com_j[, 1] %in% exclude),, drop = FALSE], 
                                                        error_val = 9999)
        }
      }
    }
    
    
    if ("FEve" %in% metrics) {
      ### Full FEve (includes all species in community) ###
      
      if (feve_calc_type == "reduced_PCoA") {
        #alternative approach to computing the minimum spanning tree that uses the same number of PCoA axes as the FRic and FDiv calculations
        diversity_data$full_FEve <- calc_feve(funcspace = x.pco_iterativedrop, com_data = com_j, mst_package = mst_package)
        
      } else {
        ### THE FD PACKAGE USES ALL AXES OF THE PCoA TO COMPUTE THE MINIMUM SPANNING TREE ###
        #to replicate the FD package method, uncomment these two lines and comment out the alternative tr.dist computation directly below
        FEve_PCoA <- x.pco_init$li[which(rownames(x.pco_init$li) %in% species_com_j),, drop = FALSE]
        #FEve_PCoA <- FEve_PCoA[match(species_com_j, rownames(FEve_PCoA)),]
        
        diversity_data$full_FEve <- calc_feve(funcspace = FEve_PCoA, com_data = com_j, mst_package = mst_package)
      }
      
      if (isTRUE(calc_contr)) {
        
        for (i in seq_len(length(species_com_j))) {
          
          ### each community constituent's contribution to FEve ###
          exclude <- species_com_j[i] #focal species to exclude from FEve calculation
          
          if (feve_calc_type == "reduced_PCoA") {
            #alternative approach to computing the minimum spanning tree that uses the same PCoA axes as the FRic and FDiv calculations
            diversity_data[i,'partial_FEve'] <- calc_feve(funcspace = x.pco_iterativedrop[!(row.names(x.pco_iterativedrop) %in% exclude),, drop = FALSE], 
                                                          com_data = com_j[ !(species_com_j %in% exclude),, drop = FALSE], 
                                                          mst_package = mst_package)
            
          } else {
            #THE APPROACH USED IN THE FD PACKAGE:
            ### THE FD PACKAGE USES ALL AXES OF THE PCoA TO COMPUTE THE MINIMUM SPANNING TREE ###
            diversity_data[i,'partial_FEve'] <- calc_feve(funcspace = FEve_PCoA[!(row.names(FEve_PCoA) %in% exclude),, drop = FALSE], 
                                                          com_data = com_j[ !(species_com_j %in% exclude),, drop = FALSE], 
                                                          mst_package = mst_package)
          }
        }
      }
    }
    
    
    ### Calculate the unstandardized and standardized functional contribution metrics ###
    if (isTRUE(calc_contr)) {
      for (m in metrics) {
        diversity_data[, paste0(m, "_prop")] <-  diversity_data[, paste0("partial_", m)]/diversity_data[, paste0("full_", m)]
        diversity_data[, paste0(m, "_contr")] <-  1 - diversity_data[, paste0(m, "_prop")]
        diversity_data[, paste0(m, "_contr_std")] <- standardize_func(diversity_data[, paste0(m, "_contr")])
      }
    }
    
    com_diversity_list[[j]] <- diversity_data #add the data/calculations for community j to com_diversity_list
    
    if (isTRUE(show_messages)) setTxtProgressBar(progress_bar, j) #update progress bar
    
  }
  
  
  #==========================================
  #== Preparation of materials for output ===
  #==========================================
  
  #collapse the dataframes within com_diversity_list into a single dataframe
  com_diversity_list_collapsed <- do.call(rbind, com_diversity_list)
  
  #add unique identifier (unique_id) if unique_id == TRUE
  if (isTRUE(calc_contr) & isTRUE(unique_id)) com_diversity_list_collapsed$unique_id <- paste(com_diversity_list_collapsed$species, com_diversity_list_collapsed$community, sep = "_")
  
  #identify the problem communities whose convex hull couldn't be formed (full_FDiv, full_FRic, partial_FDiv, and/or partial_FRic = 9999)
  included_cols <- c('full_FDiv', 'partial_FDiv', 'full_FRic', 'partial_FRic')[c('full_FDiv', 'partial_FDiv', 'full_FRic', 'partial_FRic') %in% colnames(com_diversity_list_collapsed)]
  if (length(included_cols) > 0) {
    problem_coms <- com_diversity_list_collapsed[apply(com_diversity_list_collapsed[, included_cols, drop = FALSE], 1, function(x) 9999 %in% x),]$community
  } else {
    problem_coms <- vector()
  }
  
  #remove problem communities
  com_diversity_list_FINAL <- com_diversity_list_collapsed[!(com_diversity_list_collapsed$community %in% problem_coms),, drop = FALSE]
  
  #report problem communities
  ifelse(length(problem_coms) == 0, 
         problem_com_info <- "No problem communities",
         problem_com_info <- problem_coms)
  
  #export functional diversity calculations dataframe and auxiliary information as a list
  if (isTRUE(output_PCoA)) {
    #include PCoA in output
    return(list(FD_dataframe = com_diversity_list_FINAL,
                distance_matrix_info = distance_matrix_list,
                correction = correction_message,
                communities_with_insufficient_richness = insufficient_com_info,
                absent_species_info = removed_species_info,
                failed_communities = problem_com_info,
                PCoA = x.pco_full))    
    
  } else {
    #don't include PCoA in output
    return(list(FD_dataframe = com_diversity_list_FINAL,
                distance_matrix_info = distance_matrix_list,
                correction = correction_message,
                communities_with_insufficient_richness = insufficient_com_info,
                absent_species_info = removed_species_info,
                failed_communities = problem_com_info))
  }
  
} #end of calc_fd_contr function



#####################
### MISCELLANEOUS ###
#####################

dist_wrapper <- function(df, method = c('dist', 'gowdis', 'gawdis'), arg_list = NULL) {
  ### check inputs ###
  method <- match.arg(method, several.ok = FALSE)
  if (!is.null(arg_list) & !is.list(arg_list))
    stop('arg_list needs to be a list or NULL (defaults are used).', call. = FALSE)
  
  ### process data inputs ###
  #`gawdis` seems to require categorical variables to be factors (and no characters)
  if (method == 'gawdis') {
    check_character <- sapply(df, is.character)
    if (any(check_character)) df[check_character] <- lapply(df[check_character], as.factor)
  }
  
  #input to dissimilarity function is either the trait dataframe or trait dataframe + extra arguments
  if (is.null(arg_list)) {
    arg_list_final <- list(x = df)
  } else {
    arg_list_final <- c(x = list(df), arg_list)
  }
  
  ### calculate dissimilarity matrix ###
  return(
    switch(method,
           dist = {do.call(stats::dist, arg_list_final)},
           gowdis = {do.call(FD::gowdis, arg_list_final)},
           gawdis = {do.call(gawdis::gawdis, arg_list_final)})
  )
}


listname2column <- function(list, column_name = 'element_name') {
  return(
    lapply(stats::setNames(nm = names(list)), function(NAME, data, column_name) {
      stats::setNames(cbind(NAME, data[[NAME]]), nm = c(column_name, colnames(data[[NAME]])))
    }, data = list, column_name = column_name)
  )
}


#convert pairwise distances in dataframe into a distance matrix (class dist)
df2dist_update <- function(df) {
  
  #all the elements in input df
  labels_vec <- unique(c(df[,1], df[,2]))
  
  #create square matrix with dimensions based on labels_vec length
  dist_mat <- matrix(0, length(labels_vec), length(labels_vec))
  dimnames(dist_mat) <- list(labels_vec, labels_vec)
  
  #positions of each 
  df$v1_pos <- match(df[,1], rownames(dist_mat)) #position of 
  df$v2_pos <- match(df[,2], rownames(dist_mat))
  
  #add all the pairwise differences 
  for (i in seq_len(nrow(df))) {
    dist_mat[df[i,'v1_pos'], df[i,'v2_pos']] <- df[i,3]
    dist_mat[df[i,'v2_pos'], df[i,'v1_pos']] <- df[i,3]
  }
  
  #return matrix converted to a distance matrix
  return(as.dist(dist_mat) )
}


reorder_dist <- function(dist_mat, new_order) {
  mat_convert <- as.matrix(dist_mat)
  
  return(
    as.dist(mat_convert[match(new_order, rownames(mat_convert)),
                        match(new_order, colnames(mat_convert))] )
  )
  
}





#################################
### DECOMPOSITION OF VARIANCE ###
#################################

var_decomp <- function(abundance_df, trait_df) {
  #paper describing method: https://doi.org/10.1111/j.2041-210X.2010.00071.x
  
  #STEP 1: calculate the fraction of total observations per species represented by each
  #        trait value; aka 1/sample size split by species
  trait_df1 <- trait_df %>% 
    group_by(species) %>% 
    mutate(rel_freq = 1/n()) %>% 
    ungroup()
  
  #STEP 2: - join the processed trait info with the abundance information
  #        - calculate Pindp and mean trait value for each species and add
  #          the info the merged dataframe
  combined_df <- left_join(trait_df1, abundance_df, by = 'species') %>% 
    mutate(Pindp = rel_freq*rel_abund) %>% 
    group_by(species) %>% 
    mutate(mean_trait = mean(trait)) %>% 
    ungroup()
  
  #STEP 3: calculate within species variance
  within_var <- combined_df %>% 
    group_by(species) %>%
    mutate(sst = (trait - mean_trait)^2) %>% 
    summarize(rel_abund = first(rel_abund),
              varx = mean(sst),
              .groups = 'drop') %>%
    mutate(witp = rel_abund*varx) %>% 
    pull(witp) %>% 
    sum()
  
  #STEP 4: calculate betweeen species variance
  #4a: calculate xcomp
  xcomp <- combined_df %>% 
    group_by(species) %>% 
    summarize(mean_trait = mean(trait), .groups = "drop") %>% 
    left_join(., abundance_df, by = 'species') %>% 
    mutate(xcomp = mean_trait*rel_abund) %>% 
    pull(xcomp) %>% 
    sum()
  
  #4b: calculate between species variance
  between_var <- combined_df %>% 
    group_by(species) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(sstp = ((mean_trait - xcomp)^2)*rel_abund) %>% 
    pull(sstp) %>% 
    sum()
  
  ### RETURN VARIANCE INFO ###
  return(
    data.frame(between_var = between_var,
               within_var = within_var,
               between_var_prop = between_var/(between_var + within_var),
               within_var_prop = within_var/(between_var + within_var)
    )
  )
}







FTD<-function(tdmat,weights=NULL,q=1){
  #https://github.com/ShanKothari/DecomposingFD
  
  ## contingency for one-species communities
  if(length(tdmat)==1 && tdmat==0){
    tdmat<-as.matrix(tdmat)
  }
  
  ## is the input a (symmetric) matrix or dist? if not...
  # if(!(class(tdmat) %in% c("matrix","dist"))){
  #   stop("distances must be class dist or class matrix")
  # } else if(class(tdmat)=="matrix" && !isSymmetric(unname(tdmat))){
  #   warning("trait matrix not symmetric")
  # } else if(class(tdmat)=="dist"){
  #   tdmat<-as.matrix(tdmat)
  # }
  
  ## is the input a (symmetric) matrix or dist? if not...
  if(!any(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if("matrix" %in% class(tdmat) && !isSymmetric(unname(tdmat))){
    warning("trait matrix not symmetric")
  } else if("dist" %in% class(tdmat)){
    tdmat<-as.matrix(tdmat)
  }
  
  if(!isTRUE(all.equal(sum(diag(tdmat)),0))){
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  
  if(max(tdmat)>1 || min(tdmat)<0){
    tdmat<-(tdmat-min(tdmat))/(max(tdmat)-min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  
  ## if no weights are provided, abundances are assumed equal
  if(is.null(weights)){
    nsp<-nrow(tdmat)
    weights<-rep(1/nsp,nsp)
  } else {
    nsp<-sum(weights>0)
  }
  
  if(!isTRUE(all.equal(sum(weights),1))){
    weights<-weights/sum(weights)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  
  tdmat.abund<-diag(weights) %*% tdmat %*% diag(weights)
  ## Here, because sum(weights)=1, the sum is here already adjusted by dividing by the 
  ## square of the number of species (if weights=NULL)
  ## or by multiplying by proportional abundances
  M<-sum(tdmat.abund)
  ## M equals Rao's Q in abundance-weighted calculations
  M.prime<-ifelse(nsp==1,0,M*nsp/(nsp-1))
  fij<-tdmat.abund/M
  
  ## calculating qHt
  ## fork -- if q=1, 1/(1-q) is undefined, so we use an analogue
  ## of the formula for Shannon-Weiner diversity
  ## if q!=1, we can calculate explicitly
  ## by definition, qHt=0 when all trait distances are zero
  if(isTRUE(all.equal(M,0))){
    qHt<-0
  } else if(q==1){
    fijlog<-fij*log(fij)
    fijlog[is.na(fijlog)]<-0
    qHt<-exp(-1*sum(fijlog))
  } else if(q==0){
    qHt<-sum(fij>0)
  } else {
    qHt<-sum(fij^q)^(1/(1-q))
  }
  
  ## getting qDT, qDTM, and qEt from qHt
  qDT<-(1+sqrt(1+4*qHt))/2
  qDTM<-1+qDT*M
  qEt<-qDT/nsp
  
  list(nsp=nsp,q=q,M=M,M.prime=M.prime,qHt=qHt,qEt=qEt,qDT=qDT,qDTM=qDTM)
}

## wrapper for the above function across multiple communities
## requires distance matrix w/ all species, scaled to 0-1
## community data matrix (communities x species)
## with species in same order as distance matrix
## if abund=T, values in community data matrix are treated as abundances
## otherwise, all are converted to presence-absence
## if match.names=T, the code will match species names across the
## trait distance matrix and comm data matrix and rearrange the latter

FTD.comm<-function(tdmat,spmat,q=1,abund=F,match.names=F){
  
  #https://github.com/ShanKothari/DecomposingFD
  
  ## is the input a (symmetric) matrix or dist? if not...
  if(!(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if(class(tdmat)=="matrix" && !isSymmetric(unname(tdmat))){
    warning("trait matrix not symmetric")
  } else if(class(tdmat)=="dist"){
    tdmat<-as.matrix(tdmat)
  }
  
  if(!isTRUE(all.equal(sum(diag(tdmat)),0))){
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  
  if(max(tdmat)>1 || min(tdmat)<0){
    tdmat<-(tdmat-min(tdmat))/(max(tdmat)-min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  n.comm<-nrow(spmat)
  if(match.names==T){
    sp.arr<-match(rownames(as.matrix(tdmat)),colnames(spmat))
    spmat<-spmat[,sp.arr]
  }
  
  ## apply FTD to each community in turn
  out<-apply(spmat,1,function(x) unlist(FTD(tdmat=tdmat,weights=x,q=q)))
  df.out<-data.frame(t(out))
  rownames(df.out)<-rownames(spmat)
  ## warning for zero-species communities
  if(sum(df.out$nsp==0)>0){
    warning("at least one community has no species")
  }
  
  nsp<-sum(colSums(spmat>0))
  ## this is Sw -- always an arithmetic mean, according to Evsey Kosman
  u.nsp<-mean(df.out$nsp)
  ## calculate mean richness, dispersion, evenness, FTD
  u.M<-sum(df.out$nsp*df.out$M)/sum(df.out$nsp)
  
  if(q==1){
    ## geometric mean -- limit of generalized mean as q->1
    u.qDT<-prod(df.out$qDT)^(1/n.comm)
  } else {
    ## generalized mean with m=1-q
    u.qDT<-(sum(df.out$qDT^(1-q))/n.comm)^(1/(1-q))
  }
  u.M.prime<-u.M*u.nsp/(u.nsp-1)
  
  ## calculate mean FTD and evenness
  u.qDTM<-1+u.qDT*u.M
  u.qEt<-u.qDT/u.nsp
  
  ## list more things
  list(com.FTD=df.out %>% rownames_to_column(var = 'community'),
       nsp=nsp,
       u.nsp=u.nsp,
       u.M=u.M,
       u.M.prime=u.M.prime,
       u.qEt=u.qEt,
       u.qDT=u.qDT,
       u.qDTM=u.qDTM)
}

comm.disp<-function(tdmat,com1,com2){
  #https://github.com/ShanKothari/DecomposingFD
  
  ## is the input a matrix or dist? if not...
  if(!(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if(class(tdmat)=="matrix" && !isSymmetric(unname(tdmat))){
    warning("trait matrix not symmetric")
  } else if(class(tdmat)=="dist"){
    tdmat<-as.matrix(tdmat)
  }
  
  if(isTRUE(all.equal(sum(com1),1))==F || isTRUE(all.equal(sum(com2),1))==F){
    com1<-com1/sum(com1)
    com2<-com2/sum(com2)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  
  ## mean distance between species in community 1 and 2
  mAB<-sum(diag(com1) %*% tdmat %*% diag(com2))
  ## correction so that dispersion=0 if communities have the same species
  mAA<-sum(diag(com1) %*% tdmat %*% diag(com1))
  mBB<-sum(diag(com2) %*% tdmat %*% diag(com2))
  dmAB<-mAB-(mAA+mBB)/2
  return(dmAB)
}


comm.disp.mat<-function(tdmat,spmat,abund=F,sp.weighted=FALSE){
  #https://github.com/ShanKothari/DecomposingFD
  n.comm<-nrow(spmat)
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  ## if any row doesn't sum to 1, coerce summation
  if(FALSE %in% sapply(rowSums(spmat),function(x) isTRUE(all.equal(x,1)))){
    spmat<-spmat/rowSums(spmat)
    warning("proportional abundances don't always sum to 1; summation to 1 forced")
  }
  
  disp.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) comm.disp(tdmat,com1=spmat[i,],com2=spmat[j,])))
  if(sp.weighted==T){
    nsp.comm<-rowSums(spmat>0)
    disp.mat.weight<-diag(nsp.comm) %*% disp.mat %*% diag(nsp.comm)
    return(disp.mat.weight)
  } else {
    return(disp.mat)
  }
}

M.gamma.pairwise<-function(tdmat,spmat,abund=F){
  #https://github.com/ShanKothari/DecomposingFD
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  ## if any row doesn't sum to 1, coerce summation
  if(FALSE %in% sapply(rowSums(spmat),function(x) isTRUE(all.equal(x,1)))){
    spmat<-spmat/rowSums(spmat)
    warning("proportional abundances don't always sum to 1; summation to 1 forced")
  }
  
  M.gamma<-function(tdmat,com1,com2){
    nsp1<-sum(com1>0)
    nsp2<-sum(com2>0)
    c.ind<-c(which(com1>0),which(com2>0))
    ## is this the correct abundance-weighted M.gamma?
    ## or should proportional abundances in less speciose communities count for more?
    ## alternative:
    ## c.abund<-unlist(c(com1[which(com1>0)],com2[which(com2>0)]))/2
    c.abund<-unlist(c(com1[which(com1>0)]*nsp1,com2[which(com2>0)]*nsp2))/(nsp1+nsp2)
    M.c<-sum(c.abund %*% tdmat[c.ind,c.ind] %*% c.abund)
    return(M.c)
  }
  
  n.comm<-nrow(spmat)
  M.gamma.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) M.gamma(tdmat,spmat[i,],spmat[j,])))
  return(M.gamma.mat)
}

M.beta.pairwise<-function(tdmat,spmat,abund=F,norm=F){
  nsp.comm<-rowSums(spmat>0)
  n.comm<-nrow(spmat)
  nsp.pair<-outer(1:n.comm,1:n.comm,function(i,j) nsp.comm[i]+nsp.comm[j])
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,abund=abund,sp.weighted=T)
  ## factor of 2 to include distance from A to B and B to A
  M.beta.mat<-2*disp.mat.weight/nsp.pair^2
  if(norm==T){
    M.beta.norm<-M.beta.mat/M.gamma.pairwise(tdmat,abund=abund,spmat)
    M.beta.norm[is.na(M.beta.norm)]<-0
    return(M.beta.norm)
  } else {
    return(M.beta.mat)
  }
}

FTD.beta<-function(tdmat,spmat,abund=F,q=1){
  #https://github.com/ShanKothari/DecomposingFD
  nsp<-sum(colSums(spmat)>0)
  St<-sum(spmat>0)
  n.comm<-nrow(spmat)
  
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,abund=abund,sp.weighted=T)
  M.beta<-sum(disp.mat.weight)/St^2
  M.beta.prime<-M.beta*n.comm/(n.comm-1)
  
  fAB<-disp.mat.weight/sum(disp.mat.weight)
  if(isTRUE(all.equal(sum(disp.mat.weight),0))){
    Ht.beta<-0
  } else if(q==1){
    fABlog<-fAB*log(fAB)
    fABlog[is.na(fABlog)]<-0
    Ht.beta<-exp(-1*sum(fABlog))
  } else if(q==0){
    Ht.beta<-sum(fAB>0)
  } else {
    Ht.beta<-sum(fAB^q)^(1/(1-q))
  }
  qDT.beta<-(1+sqrt(1+4*Ht.beta))/2
  qDTM.beta<-1+qDT.beta*M.beta
  Et.beta<-qDT.beta/n.comm
  
  list(nsp=nsp,St=St,n.comm=n.comm,q=q,M.beta=M.beta,M.beta.prime=M.beta.prime,Ht.beta=Ht.beta,qDT.beta=qDT.beta,qDTM.beta=qDTM.beta,disp.mat.weight=disp.mat.weight)
}