#### random intercep BART model fitting function
riBMfits <- function(data, # A matrix or data frame possibly containing confounder(s), exposure(s), outcome(s) and mortality indicator(s) in temporal order from left to right.
                     var.type, # Vector of variable specifications for data. X=confounder, Fi=fixed (e.g. the exposure), "Ra"=random exposure, Y=outcome, S=survival, D= dropout indicator.
                     linear_formula, ##### A model formula with the linear variables on the right-hand-side (left-hand-side is not used).
                     fixed.regime = NULL, # A list specifying the possibly contrasting regime(s) for exposed/treated and unexposed/untreated. Can also be used for fixing a confounder.
                     random.regime = NULL, # "uniform", "normal" or "skew.normal" default is NULL
                     param = NULL, # a vector of length 1 (in the binomial case), 2 (in the uniform/normal case) or 3 (in the skew normal case),
                     # whose components represent the parameters linked to a binomial(a), uniform(a, b), a normal(a, b) or a skew normal(a, b, c)
                     above = FALSE,
                     cutoff = NULL,
                     nat_value = FALSE,
                     incremental = FALSE,
                     drop_param = NULL,
                     J = 2000, # Size of pseudo data. Default is set to 2,000.
                     opts = NULL, # opts see SoftBart
                     Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                     By = 100, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                     weighted = FALSE, # if TRUE overlap weights are computed and the a weighted outcome is returned
                     tgroup = rep(1,  ncol(data) - 1),
                     num_tree = 20,
                     alpha = 1,
                     eta = 1,
                     phi = rep(1, max(tgroup) - 1),
                     alpha_vec = rep(1, max(tgroup) - 1),
                     alpha_shape_1 = 0.5,
                     remove_intercept = TRUE,
                     print_iter = FALSE,
                     ...
                     
) {
  
  #### the covariates in linear function.
  Z <- model.matrix(linear_formula, data)
  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }
  data <- data[,!(names(data) %in% colnames(Z) ) ] ## exclude covariates from data
  
  data <- as.matrix(data)
  
  
  ## test if the variables in data are continuous or binary etc
  continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1))
  
  # tmp <- numeric(ncol(data))
  # for(i in 1:ncol(data)) {
  #     if(i == 1) {
  #         tmp[1] <- all(na.omit(data[,1]) == 0) | all(na.omit(data[,1]) == 1)
  #     } else {
  #         tmp_dat <- na.omit(data[,1:i])
  #         tmp[i] <- all(tmp_dat[,i] == 0) | all(tmp_dat[,i] == 1)
  #     }
  # }
  #
  # data <- data[,tmp==0]
  # var.type <- var.type[tmp==0]
  # tgroup <- tgroup[tmp[-length(tmp)]==0]
  
  n_Y <- length(which(var.type=="Y")) # Number of outcome variables, >1 allowed
  n_S <- length(which(var.type=="S")) # Number of survival variables, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  n_Ra <- length(which(var.type=="Ra"))
  
  
  if(!is.null(fixed.regime)) {
    n_Reg <- max(ifelse(!is.null(fixed.regime), length(fixed.regime), 0))
  } else {
    n_Reg <- 1
  }
  
  if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  if (any(is.na(data[, which(var.type == "X0")]))) stop("The baseline covariates includes NAs ")
  
  if(!is.null(fixed.regime)){
    if (all(sapply(fixed.regime, length) != n_Fi)) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  ##################################################################3
  ## Fit the observed data models using BART and save as a list
  BModels <- vector("list", ncol(data) - 1)
  if(weighted == TRUE) {
    WModels <- vector("list", n_Ra)
  }
  if(is.null(opts)){
    opts <- Opts(num_burn = num_save, num_thin = num_thin, num_save = num_save,
                 update_sigma_mu = TRUE, update_s = TRUE, update_alpha = TRUE, update_beta = FALSE,
                 update_gamma = FALSE, update_tau = TRUE, update_tau_mean = FALSE,
                 update_sigma = TRUE, cache_trees = TRUE)
  }
  we <- 1
  for (i in 1:(ncol(data)-1)) {
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
    if(i==1){
      if(var.type[1+i] == "Fi" | var.type[1+i] == "X0" ) {
        BModels[[i]]  <- NA
      } else if(var.type[1+i] == "Ra") {
        if(all(nat_value == FALSE)){
          BModels[[i]]  <- NA
        } else {
          X_train <- as.matrix(cbind(1, data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D"), drop=FALSE]))
          id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
          Z_train <- Z[id == 1, id_col == 1]
          Y_train <- data[id == 1, (1 + i)]
          
          ## group static covariates w/ last time point
          tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
          tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
          
          #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
          if(abs(sum(Y_train, na.rm = TRUE)) < 2){
            new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
            hypers <- Hypers(
              X = cbind(X_train, Z_train), Y = Y_train,
              tgroup = tmp_tg, num_tree = num_tree,
              alpha = alpha, eta = eta, phi = phi,
              alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
              sigma_hat = new_sig,            # <-- bypass glmnet
              normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
            ) }else{
              hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                               tgroup = tmp_tg,
                               num_tree = num_tree,
                               alpha = alpha,
                               eta = eta, phi = phi,
                               alpha_vec = alpha_vec,
                               alpha_shape_1 = alpha_shape_1) }
          
          if (Suppress == TRUE) {
            quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
          } else if (Suppress == FALSE) {
            BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
          }
        }
      } else {
        X_train <- as.matrix(cbind(1, data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D"), drop=FALSE]))
        id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
        Z_train <- Z[id == 1, id_col == 1]
        Y_train <- data[id == 1, (1 + i)]
        tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
        tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
        #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
        if(abs(sum(Y_train, na.rm = TRUE)) < 2){
          new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
          hypers <- Hypers(
            X = cbind(X_train, Z_train), Y = Y_train,
            tgroup = tmp_tg, num_tree = num_tree,
            alpha = alpha, eta = eta, phi = phi,
            alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
            sigma_hat = new_sig,            # <-- bypass glmnet
            normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
          ) }else{
            hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                             tgroup = tmp_tg,
                             num_tree = num_tree,
                             alpha = alpha,
                             eta = eta, phi = phi,
                             alpha_vec = alpha_vec,
                             alpha_shape_1 = alpha_shape_1) }
        if (continuous[1+i] == FALSE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == TRUE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == FALSE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers)
        } else if (continuous[1+i] == TRUE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
        }
      }
    } else if(i==2){
      if(var.type[1+i] == "Fi" | var.type[1+i] == "X0" ) {
        BModels[[i]]  <- NA
      } else if(var.type[1+i] == "Ra") {
        if(nat_value == FALSE) {
          BModels[[i]]  <- NA
        } else {
          X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
          id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
          Z_train <- Z[id == 1, id_col == 1]
          Y_train <- data[id == 1, (1 + i)]
          tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
          tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
          #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
          if(abs(sum(Y_train, na.rm = TRUE)) < 2){
            new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
            hypers <- Hypers(
              X = cbind(X_train, Z_train), Y = Y_train,
              tgroup = tmp_tg, num_tree = num_tree,
              alpha = alpha, eta = eta, phi = phi,
              alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
              sigma_hat = new_sig,            # <-- bypass glmnet
              normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
            ) }else{
              hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                               tgroup = tmp_tg,
                               num_tree = num_tree,
                               alpha = alpha,
                               eta = eta, phi = phi,
                               alpha_vec = alpha_vec,
                               alpha_shape_1 = alpha_shape_1) }
          if (Suppress == TRUE) {
            quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
          } else if (Suppress == FALSE) {
            BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
          }
        }
        
        if(weighted == TRUE) {
          ## compute overlap-weights models
          if(above == TRUE){
            y_tmp <- ifelse(data[id == 1, (1 + i)] < cutoff[we], 1, 0)
          } else {
            y_tmp <- ifelse(data[id == 1, (1 + i)] >  cutoff[we], 1, 0)
          }
          if(all(y_tmp == 1) | all(y_tmp == 0)) {
            WModels[[we]]  <- NULL
          } else {
            X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
            id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
            Z_train <- Z[id == 1, id_col == 1]
            Y_train <- y_tmp
            tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
            tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
            #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
            if(abs(sum(Y_train, na.rm = TRUE)) < 2){
              new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
              hypers <- Hypers(
                X = cbind(X_train, Z_train), Y = Y_train,
                tgroup = tmp_tg, num_tree = num_tree,
                alpha = alpha, eta = eta, phi = phi,
                alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
                sigma_hat = new_sig,            # <-- bypass glmnet
                normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
              ) }else{
                hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                                 tgroup = tmp_tg,
                                 num_tree = num_tree,
                                 alpha = alpha,
                                 eta = eta, phi = phi,
                                 alpha_vec = alpha_vec,
                                 alpha_shape_1 = alpha_shape_1) }
            if (Suppress == TRUE) {
              quiet(WModels[[we]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers))
            } else  if (Suppress == FALSE) {
              WModels[[we]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers)
            }
          }
          we <- we + 1
        }
      } else {
        X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
        id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
        Z_train <- Z[id == 1, id_col == 1]
        Y_train <- data[id == 1, (1 + i)]
        tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
        tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
        #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
        if(abs(sum(Y_train, na.rm = TRUE)) < 2){
          new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
          hypers <- Hypers(
            X = cbind(X_train, Z_train), Y = Y_train,
            tgroup = tmp_tg, num_tree = num_tree,
            alpha = alpha, eta = eta, phi = phi,
            alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
            sigma_hat = new_sig,            # <-- bypass glmnet
            normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
          ) }else{
            hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                             tgroup = tmp_tg,
                             num_tree = num_tree,
                             alpha = alpha,
                             eta = eta, phi = phi,
                             alpha_vec = alpha_vec,
                             alpha_shape_1 = alpha_shape_1) }
        if (continuous[1+i] == FALSE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == TRUE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == FALSE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers)
        } else if (continuous[1+i] == TRUE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
        }
      }
    } else if(i>2) {
      if(var.type[1+i] == "Fi" | var.type[1+i] == "X0" ) {
        BModels[[i]]  <- NA
      } else if(var.type[1+i] == "Ra") {
        if(nat_value == FALSE) {
          BModels[[i]]  <- NA
        } else {
          X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
          id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
          Z_train <- Z[id == 1, id_col == 1]
          Y_train <- data[id == 1, (1 + i)]
          tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
          tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
          #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
          if(abs(sum(Y_train, na.rm = TRUE)) < 2){
            new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
            hypers <- Hypers(
              X = cbind(X_train, Z_train), Y = Y_train,
              tgroup = tmp_tg, num_tree = num_tree,
              alpha = alpha, eta = eta, phi = phi,
              alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
              sigma_hat = new_sig,            # <-- bypass glmnet
              normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
            ) }else{
              hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                               tgroup = tmp_tg,
                               num_tree = num_tree,
                               alpha = alpha,
                               eta = eta, phi = phi,
                               alpha_vec = alpha_vec,
                               alpha_shape_1 = alpha_shape_1) }
          if (Suppress == TRUE) {
            quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
          } else if (Suppress == FALSE) {
            BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
          }
        }
        
        if(weighted == TRUE) {
          ## compute overlap-weights models
          if(above == TRUE){
            y_tmp <- ifelse(data[id == 1, (1 + i)] < cutoff[we], 1, 0)
          } else {
            y_tmp <- ifelse(data[id == 1, (1 + i)] >  cutoff[we], 1, 0)
          }
          if(all(y_tmp == 1) | all(y_tmp == 0)) {
            WModels[[we]]  <- NULL
          } else {
            X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
            id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
            Z_train <- Z[id == 1, id_col == 1]
            Y_train <- y_tmp
            tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
            tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
            #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
            if(abs(sum(Y_train, na.rm = TRUE)) < 2){
              new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
              hypers <- Hypers(
                X = cbind(X_train, Z_train), Y = Y_train,
                tgroup = tmp_tg, num_tree = num_tree,
                alpha = alpha, eta = eta, phi = phi,
                alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
                sigma_hat = new_sig,            # <-- bypass glmnet
                normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
              ) }else{
                hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                                 tgroup = tmp_tg,
                                 num_tree = num_tree,
                                 alpha = alpha,
                                 eta = eta, phi = phi,
                                 alpha_vec = alpha_vec,
                                 alpha_shape_1 = alpha_shape_1) }
            if (Suppress == TRUE) {
              quiet(WModels[[we]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train,  opts = opts, hypers = hypers))
            } else  if (Suppress == FALSE) {
              WModels[[we]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers)
            }
          }
          we <- we + 1
        }
      } else {
        X_train <- data[id == 1, which(var.type[1:i]!="S" & var.type[1:i]!="D")]
        id_col <- ifelse(apply(Z[id == 1, ], 2, function(x) (!all(x == 0) )), 1, 0) ### select not all 0 in countries
        Z_train <- Z[id == 1, id_col == 1]
        Y_train <- data[id == 1, (1 + i)]
        tmp_tg <- tgroup[which(var.type[1:i]!="S" & var.type[1:i]!="D")]
        tmp_tg[which(tmp_tg==0)] <- max(tmp_tg)
        #### There is one position of some Y(s) is 1, others are 0, that the glmnet() doesn't work when applied Hypers(), hereby we give a value of sigma_aht 
        if(abs(sum(Y_train, na.rm = TRUE)) < 2){
          new_sig <- sqrt(mean(Y_train) * (1 - mean(Y_train)) + 1e-8)      # stable default
          hypers <- Hypers(
            X = cbind(X_train, Z_train), Y = Y_train,
            tgroup = tmp_tg, num_tree = num_tree,
            alpha = alpha, eta = eta, phi = phi,
            alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
            sigma_hat = new_sig,            # <-- bypass glmnet
            normalize_Y = FALSE         # avoid extra scaling for extreme imbalance
          ) }else{
            hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train,
                             tgroup = tmp_tg,
                             num_tree = num_tree,
                             alpha = alpha,
                             eta = eta, phi = phi,
                             alpha_vec = alpha_vec,
                             alpha_shape_1 = alpha_shape_1) }
        if (continuous[1+i] == FALSE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == TRUE & Suppress == TRUE) {
          quiet(BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers))
        } else if (continuous[1+i] == FALSE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_probit(X = X_train, Y = Y_train,  Z = Z_train, opts = opts, hypers = hypers)
        } else if (continuous[1+i] == TRUE & Suppress == FALSE) {
          BModels[[i]]  <- gc_gsoftbart_regression(X = X_train, Y = Y_train, Z = Z_train, opts = opts, hypers = hypers)
        }
      }
    }
    if(print_iter == TRUE){print(i)}
  }
  return(BModels)
}

#### G-computation for the random intercep softbart model

gcomp_risoftbart <- function(data, # A matrix or data frame possibly containing confounder(s), exposure(s), outcome(s) and mortality indicator(s) in temporal order from left to right.
                             var.type, # Vector of variable specifications for data. X=confounder, Fi=fixed (e.g. the exposure), "Ra"=random exposure, Y=outcome, S=survival, D= dropout indicator.
                             linear_formula, ##### A model formula with the linear variables on the right-hand-side (left-hand-side is not used).
                             fixed.regime = NULL, # A list specifying the possibly contrasting regime(s) for exposed/treated and unexposed/untreated. Can also be used for fixing a confounder.
                             random.regime = NULL, # "uniform", "normal" or "skew.normal" default is NULL
                             param = NULL, # a vector of length 1 (in the binomial case), 2 (in the uniform/normal case) or 3 (in the skew normal case),
                             # whose components represent the parameters linked to a binomial(a), uniform(a, b), a normal(a, b) or a skew normal(a, b, c)
                             above = FALSE,
                             cutoff = NULL,
                             nat_value = FALSE,
                             incremental = FALSE,
                             drop_param = NULL,
                             J = 2000, # Size of pseudo data. Default is set to 2,000.
                             opts = NULL, # opts see SoftBart
                             Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                             By = 100, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                             weighted = FALSE, # if TRUE overlap weights are computed and the a weighted outcome is returned
                             tgroup = rep(1,  ncol(data) - 1),
                             num_tree = 20,
                             alpha = 1,
                             eta = 1,
                             phi = rep(1, max(tgroup) - 1),
                             alpha_vec = rep(1, max(tgroup) - 1),
                             alpha_shape_1 = 0.5,
                             BModels = NULL,
                             time.type = NULL, #### A vector describes the time of variables from the argument var.type, this only use to predict the missing random intercept(s).
                             impute_method = NULL, # BART or MVN
                             ridge = 1e-6,
                             remove_intercept = TRUE,
                             ...
                             
) 
{
  
  #### the covariates in linear function.
  Z <- model.matrix(linear_formula, data)
  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }
  data <- data[,!(names(data) %in% colnames(Z) ) ] ## exclude covariates from data
  
  data <- as.matrix(data)
  
  
  ## test if the variables in data are continuous or binary etc
  continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1))
  
  # tmp <- numeric(ncol(data))
  # for(i in 1:ncol(data)) {
  #     if(i == 1) {
  #         tmp[1] <- all(na.omit(data[,1]) == 0) | all(na.omit(data[,1]) == 1)
  #     } else {
  #         tmp_dat <- na.omit(data[,1:i])
  #         tmp[i] <- all(tmp_dat[,i] == 0) | all(tmp_dat[,i] == 1)
  #     }
  # }
  #
  # data <- data[,tmp==0]
  # var.type <- var.type[tmp==0]
  # tgroup <- tgroup[tmp[-length(tmp)]==0]
  
  
  n_Y <- length(which(var.type=="Y")) # Number of outcome variables, >1 allowed
  n_S <- length(which(var.type=="S")) # Number of survival variables, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  n_Ra <- length(which(var.type=="Ra"))
  
  if(!is.null(fixed.regime)) {
    n_Reg <- max(ifelse(!is.null(fixed.regime), length(fixed.regime), 0))
  } else {
    n_Reg <- 1
  }
  
  if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  if (any(is.na(data[, which(var.type == "X0")]))) stop("The baseline covariates includes NAs ")
  
  if(!is.null(fixed.regime)){
    if (all(sapply(fixed.regime, length) != n_Fi)) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  ##################################################################3
  ## Fit the observed data models using General Soft BART and save as a list
  if(is.null(BModels)){
    stop("BM shouldn't be NULL, which implies there is no fitted model for prediction.")
  } 
  ########################################################
  ## Create matrices to store the output. One for the outcome and one for survival probabilities.
  ## Note, the code allows for multiple outcomes/survival, e.g. for different time points.
  
  y <- matrix(nrow = ifelse(!is.null(fixed.regime) | !is.null(random.regime), n_Reg, 1) * n_Y, ncol = opts$num_save)
  s <- matrix(nrow = ifelse(!is.null(fixed.regime) | !is.null(random.regime), n_Reg, 1) * n_S, ncol = opts$num_save)
  if(weighted == TRUE) {
    yw <- matrix(nrow = ifelse(!is.null(fixed.regime) | !is.null(random.regime), n_Reg, 1) * n_Y, ncol = opts$num_save)
  } else {
    yw <- NULL
  }
  if(weighted == TRUE) {
    sw <- matrix(nrow = ifelse(!is.null(fixed.regime) | !is.null(random.regime), n_Reg, 1) * n_S, ncol = opts$num_save)
  } else {
    sw <- NULL
  }
  
  iterations <- seq(from = opts$num_burn + opts$num_thin,
                    to = opts$num_burn + opts$num_thin * opts$num_save,
                    by = opts$num_thin)
  
  ## MC-integration over "opts$num_save" iterations
  for(it in 1:opts$num_save) {
    ks <- 1
    ky <- 1
    l <- 1
    n <- 1
    m <- 1
    o <- 1
    w <- 1
    d <- 0
    drop <- NULL
    s_hat <- NULL
    w_hat <- NULL
    # sample data for the first confounder
    #if(var.type[1] == "Fi") {
    #  x <- Map(matrix, lapply(fixed.regime, function(r) rep(r[1], J)))
    #  l <- l + 1
    #  m <- 2
    #  } else
    
    new_BModels <- BModels #### To save the imputed groups effects for each iteration, but not affects other iterations.
    
    if(var.type[1] == "Ra") {
      if(weighted == TRUE){
        if(is.null(WModels[[1]])){
          w_hat <- 1
        } else {
          #### ???? can not find the x before this line
          w_hat <- lapply(x, function(x2) { ri_pred_w(WModels[[1]], x2, iterations[it],
                                                   remove_intercept) } )
        }
        w <- w + 1
      }
      if(random.regime[1] == "uniform") {
        x <- list(matrix(runif(J, param[[1]][1], param[[1]][2])))
      } else if (random.regime[1] == "normal") {
        x <- list(matrix(rnorm(J, param[[1]][1], param[[1]][2])))
      } else if(random.regime[1] == "triangular"){
        x <- list(matrix(EnvStats::rtri(J, param[[1]][1], param[[1]][2], param[[1]][3])))
      } else if (random.regime[1] == "binomial") {
        x <- list(matrix(rbinom(J, 1, param[[1]][1])))
      }
      m <- 2
    } else if(var.type[1] == "X0") {
      dataC0 <- data.frame(cbind(data[, which(var.type == "X0")], Z) )
      x <- list(BayesBoot(dataC0, J))
      m <- ncol(dataC0) - ncol(Z) + 1
    } else if(var.type[1] == "Fi" & var.type[2] == "X0") {
      data0 <- Map(data.frame, lapply(fixed.regime, function(r) cbind(data[which(data[,1]==r[1]), c(1, which(var.type == "X0"))], Z) ))
      x <- lapply(data0, function(x) BayesBoot(x, J))
      m <- ncol(data0[[1]]) - ncol(Z) + 1
    } else {
      dataC0 <- data.frame(cbind(data[, 1, drop = FALSE], Z))
      x <- list(BayesBoot(dataC0, J))
      m <- 2
    }
    
    for (j in m:ncol(data) ) { 
      
      if(is.null(time.type)) {
        current_BModels <- new_BModels[[j - 1]]
      } else {
        
        ### Identify the time of column i 
        current_time <- as.integer(sub("^T", "", as.character(time.type)))[j] 
        if(var.type[j] == "Fi"){
          current_BModels <- new_BModels[[j - 1]]
        }else { ### For "S", "D", "X", "Y"
          
          ### Selection of the indices of prior models
          prev_idx <- select_history_indices(j, var.type, time.type) 
          if(length(prev_idx) == 0) {
            current_BModels <- pred_beta_BM_noHis(new_BModels[[j - 1]], linear_formula, iterations[it])
            new_BModels[[j - 1]] <- current_BModels ### update the model
          }else{ 
             ### For "X", previous models might contain the model of "X0" 
            if(any(var.type[prev_idx] == "X0") == TRUE & length(prev_idx) == 1){
              current_BModels <- pred_beta_BM_noHis(new_BModels[[j - 1]], linear_formula, iterations[it])
              new_BModels[[j - 1]] <- current_BModels ### update the model
              
            }else if(any(var.type[prev_idx] == "X0") == TRUE & length(prev_idx) > 1){
             
               new_prev_idx <- prev_idx[var.type[prev_idx] != "X0"]
              
              prev_BModels <- list()
              for (n_idx in seq_along(new_prev_idx)) {
                prev_BModels[[n_idx]] <- new_BModels[[new_prev_idx[n_idx]-1]]
              }
              
                current_BModels <- pred_beta_BM(new_BModels[[j - 1]], 
                                                prev_BModels, 
                                                linear_formula, 
                                                iterations[it],
                                                impute_method = impute_method,
                                                ridge)
          
              
              new_BModels[[j - 1]] <- current_BModels ### update the model
            }else{  ##  for all "X" history
              
              prev_BModels <- list()
              for (n_idx in seq_along(prev_idx)) {
                prev_BModels[[n_idx]] <- new_BModels[[prev_idx[n_idx]-1]]
              }
            #current_BModels <- pred_beta_BM(new_BModels[[j - 1]], prev_BModels, linear_formula, iterations[it])
              current_BModels <- pred_beta_BM(new_BModels[[j - 1]], 
                                              prev_BModels, 
                                              linear_formula, 
                                              iterations[it],
                                              impute_method = impute_method,
                                              ridge)
              
              new_BModels[[j - 1]] <- current_BModels ### update the model
            }
          }  
        }
      } 
      
      
      
      
      if(var.type[j] == "X") {
        ##### check which the random intercept parameters of the specific groups miss at j
        x <- lapply(x, function(x2) {ri_pred_comb(continuous[j], current_BModels, x2, iterations[it], remove_intercept)})
      } else if(var.type[j] == "Fi") {
        x <- Map(cbind, x, lapply(fixed.regime, function(r) r[l]))
        l <- l + 1
      } else if(var.type[j] == "Ra") {
        if(weighted == TRUE) {
          if(is.null(WModels[[w]])){
            w_hat_tmp <- 1
          } else {
            w_hat_tmp <- lapply(x, function(x2) {ri_pred_w(WModels[[w]], x2, iterations[it], remove_intercept) } )
          }
          if(!is.null(w_hat)){
            w_hat <- Map(function(a, b){a*b}, w_hat, w_hat_tmp)
          } else {
            w_hat <- w_hat_tmp
          }
          w <- w + 1
        }
        x <- lapply(x, function(x2) { ri_pred_ra(continuous[j], current_BModels,
                                              x2, iterations[it], random.regime[o], param[[o]],
                                              above, cutoff[o], nat_value, incremental,
                                              remove_intercept) } )
        o <- o + 1
      } else if(var.type[j] == "S") {
        surv_hat <- lapply(x, function(x2) {ri_pred_surv(current_BModels, x2, iterations[it],
                                                      remove_intercept)})
        if(!is.null(s_hat)){
          s_tmp <- Map(function(a, b){a*b}, s_hat, surv_hat)
        } else {
          s_tmp <- surv_hat
        }
        s[n:((n-1) + n_Reg), it] <- sapply(s_tmp, mean)
        n <- n + n_Reg
        surv <- lapply(surv_hat, function(a){rbinom(length(a), 1, prob = a)})
        x <- Map(function(a, b){a[b == 1,]}, x, surv)
        s_hat <- Map(function(a, b){a[b == 1]}, s_tmp, surv)
        if(!is.null(w_hat)){
          sw_hat <- Map(function(a, b){(a * b)/sum(b)}, s_tmp, w_hat)
          sw[ks:((ks-1) + n_Reg), it] <- sapply(sw_hat, sum)
          w_hat <- Map(function(a, b){a[b == 1]}, w_hat, surv)
        } else if(is.null(w_hat) & weighted == TRUE){
          sw[ks:((ks-1) + n_Reg), it] <- sapply(s_tmp, mean)
        }
        ks <- ks + n_Reg
        if(!is.null(drop)){
          drop <- Map(function(a, b){a[b == 1]}, drop, surv)
        }
      } else if(var.type[j] == "D") {
        drop_tmp <- lapply(x, function(x2) {ri_pred_drop(current_BModels, x2, iterations[it], remove_intercept)})
        if(d == 0){
          drop <- drop_tmp
        } else {
          drop <- Map(function(e, f){ifelse(e==0 & f==0, 0, 1)}, drop_tmp, drop)
        }
        d <- d + 1
        if(is.null(drop_param)){
          drop_shift <- lapply(drop, function(f){f*0})
        } else {
          drop_shift <- lapply(drop, function(f){f * EnvStats::rtri(length(f), drop_param[[d]][1], drop_param[[d]][2], drop_param[[d]][3])})
        }
      } else if(
        var.type[j] == "Y") {
        tmp <- lapply(x, function(x2) {
          ri_pred_y(continuous[j], current_BModels, x2, iterations[it], remove_intercept)} )
        
        if(d > 0){
          tmp <- Map(function(e, f){e + f}, tmp, drop_shift)
        }
        y[ky:((ky-1) + n_Reg), it] <- sapply(tmp, mean)
        if(!is.null(w_hat)){
          yw_hat <- Map(function(a, b){(a * b)/sum(b)}, tmp, w_hat)
          yw[ky:((ky-1) + n_Reg), it] <- sapply(yw_hat, sum)
        }
        ky <- ky + n_Reg
        if(continuous[j] == FALSE) {
          tmp <- lapply(tmp, function(tmp2) rbinom(length(tmp2), 1, prob = tmp2))
          x <- Map(cbind, x, tmp)
        } else if(continuous[j] == TRUE) {
          tmp <- lapply(tmp, function(tmp2) rnorm(length(tmp2), mean = tmp2, sd = mean(current_BModels$sigma)))
          x <- Map(cbind, x, tmp)
        }
      }
    }
    if(it %in% seq(0, opts$num_save, by = By)) {
      print(paste("done ", it, " (out of ", opts$num_save,")", sep=""))
    }
  }
  
  
  ## save the posterior samples of the mean predicted the outcome(s) (y_hat) and summaries
  if(n_Y > 0){
    means_out <- t(apply(y, 1, function(x) {c(mean(x), quantile(x, c(0.025, 0.975), na.rm=TRUE))}))
    if(is.null(fixed.regime)){
      rownames(means_out) <- paste0("Y ", 1:n_Y)
      rownames(y) <- paste0("Y ", 1:n_Y)
      
    } else {
      rownames(means_out) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y ", rep(1:n_Y, each=n_Reg))
      rownames(y) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y ", rep(1:n_Y, each=n_Reg))
    }
    
    colnames(means_out) <- c("mu", "2.5%", "97.5%")
  }  else {
    y <- NA
    means_out <- NA
  }
  if(n_Y > 0 & weighted == TRUE){
    means_wout <- t(apply(yw, 1, function(x) {c(mean(x), quantile(x, c(0.025, 0.975), na.rm=TRUE))}))
    if(is.null(fixed.regime)){
      rownames(means_out) <- paste0("Y ", 1:n_Y)
      rownames(y) <- paste0("Y ", 1:n_Y)
      
    } else {
      rownames(means_wout) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y, weighted ", rep(1:n_Y, each=n_Reg))
      rownames(yw) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y, weighted ", rep(1:n_Y, each=n_Reg))
    }
    
    colnames(means_wout) <- c("mu", "2.5%", "97.5%")
  }  else {
    yw <- NA
    means_wout <- NA
  }
  ## save the posterior samples of the mean survival probabilities and summaries
  if(n_S>0){
    means_surv <- t(apply(s, 1, function(x) {c(mean(x), quantile(x, c(0.025, 0.975), na.rm=TRUE))}))
    if(is.null(fixed.regime)){
      rownames(means_surv) <- paste0("Survival ", 1:n_S)
      rownames(s) <- paste0("Survival ", 1:n_S)
      
    } else {
      rownames(means_surv) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival ", rep(1:n_S, each=n_Reg))
      rownames(s) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival ", rep(1:n_S, each=n_Reg))
    }
    colnames(means_surv) <- c("mu", "2.5%", "97.5%")
  } else {
    s <- NA
    means_surv <- NA
  }
  if(n_S>0 & weighted == TRUE){
    means_wsurv <- t(apply(sw, 1, function(x) {c(mean(x), quantile(x, c(0.025, 0.975), na.rm=TRUE))}))
    if(is.null(fixed.regime)){
      rownames(means_wsurv) <- paste0("Survival, weighted ", 1:n_S)
      rownames(sw) <- paste0("Survival, weighted ", 1:n_S)
      
    } else {
      rownames(means_wsurv) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival, weighted ", rep(1:n_S, each=n_Reg))
      rownames(sw) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival, weighted ", rep(1:n_S, each=n_Reg))
    }
    colnames(means_wsurv) <- c("mu", "2.5%", "97.5%")
  } else {
    sw <- NA
    means_wsurv <- NA
  }
  
  
  ####
  out <- list(means_out, y, means_wout, yw, means_surv, s, means_wsurv, sw)
  names(out) <- list("summary_out", "y_hat", "summary_weighted_out",
                     "weighted_y_hat", "summary_surv", "s_hat",
                     "summary_weighted_surv", "weighted_s_hat")
  
  return(out)
}

