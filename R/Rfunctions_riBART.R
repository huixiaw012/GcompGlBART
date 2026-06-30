BayesBoot <- function(df , J) {
  BBweights <- rBeta2009::rdirichlet(1, rep(1, nrow(df)))
  BBdraws <- rmultinom(1, J, BBweights)
  df0 <- as.matrix(data.frame(lapply(df, rep, BBdraws))) #repeat i^th row of input_df, BBdraw[i] times
  return(df0)
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

ri_pred_comb <- function(cont, BM, MCdatatmp, IT, remove_intercept){

  if(ncol(MCdatatmp)==1) {
    MCdata <- as.data.frame(cbind(1, MCdatatmp))
  } else {
    MCdata <- as.data.frame(MCdatatmp)
  }

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula(
    paste("~", paste(present_vars, collapse = " + "))
  ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]


  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(opts$num_save == 1){
    beta <- BM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(BM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]

  if(cont == FALSE) {
    #x_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset) ### need to change
    #new_MCdata <- cbind(MCdata, rbinom(length(x_hat), 1, prob = x_hat)) ### need to change
    r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) + BM$offset
    eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
    x_z_hat <- r + eta
    new_MCdata <- cbind(MCdata, rbinom(length(x_z_hat), 1, prob =pnorm(x_z_hat)))

  } else if(cont == TRUE) {

    r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) * BM$sd_Y + BM$mu_Y
    eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
    x_z_hat <- r + eta
    #x_hat <- as.numeric(BM$forest$predict_iteration(X, IT)) * BM$sd_Y + BM$mu_Y
    new_MCdata <- cbind(MCdata, rnorm(length(x_z_hat), mean = x_z_hat, sd = mean(BM$sigma)))
    #cbind(X, rnorm(length(x_z_hat), mean = x_z_hat, sd = mean(BM$sigma)), Z)

  }
  return(new_MCdata)
}

ri_pred_drop <- function(BM, MCdatatmp, IT, remove_intercept){

  #### full_linear_formula : the linear formula of new dataset at baseline.
  #### ln.formula : the pratical exited linear formula of fitting model.
  MCdata <- as.data.frame(MCdatatmp)

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula(
    paste("~", paste(present_vars, collapse = " + "))
  ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]


  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(opts$num_save == 1){
    beta <- BM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(BM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]

  r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) + BM$offset
  eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
  x_z_hat <- r + eta
  dropout_ind <- rbinom(length(x_z_hat), 1, prob =pnorm(x_z_hat))
  #x_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
  #dropout_ind <- rbinom(length(x_hat), 1, prob = x_hat)

  return(dropout_ind)
}


ri_pred_ra <- function(cont, BM, MCdatatmp, IT, random.regime, param, above, cutoff, nat_value, incremental, remove_intercept){ #, rank){ # above = NULL
  if(ncol(MCdatatmp)==1) {
    MCdata <- as.data.frame(cbind(1, MCdatatmp))
  } else {
    MCdata <- as.data.frame(MCdatatmp)
  }

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula(
    paste("~", paste(present_vars, collapse = " + "))
  ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]


  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(opts$num_save == 1){
    beta <- BM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(BM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]


  if(nat_value == TRUE){ # based on the natural value
    if(cont == FALSE) {
      r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) + BM$offset
      eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
      fi_hat <- pnorm(r + eta)
      #fi_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
      new_fi <- rbinom(length(fi_hat), 1, prob = fi_hat)
    } else if(cont == TRUE) {
      r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) * BM$sd_Y + BM$mu_Y
      eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
      fi_hat <- pnorm(r + eta)
      #fi_hat <- as.numeric(BM$forest$predict_iteration(X, IT)) * BM$sd_Y + BM$mu_Y
      new_fi <- rnorm(length(fi_hat), mean = fi_hat, sd = mean(BM$sigma))
      if(above == TRUE){
        interv <- which(new_fi > cutoff)
      } else if (above == FALSE) {
        interv <- which(new_fi < cutoff)
      }
      if(random.regime == "uniform") {
        int_dist <- runif(length(interv), param[1], param[2])
      } else if (random.regime == "normal") {
        int_dist <- rnorm(length(interv), param[1], param[2])
      } else if(random.regime == "triangular"){
        int_dist <- EnvStats::rtri(length(interv), param[1], param[2], param[3])
      } else if (random.regime == "binomial") {
        int_dist <- rbinom(length(interv), 1, param[1])
      }
      if (incremental == TRUE){
        new_fi[interv] <- new_fi[interv] + int_dist
      } else {
        new_fi[interv] <- int_dist
      }

      #if(rank == TRUE){
      #  new_fi[interv][order(new_fi[interv])] <- int_dist[order(int_dist)]
      #} else {
      #new_fi[interv] <- int_dist
      #}
    }
  } else if(nat_value == FALSE){ # not based on the natural value
    if(random.regime == "uniform") {
      int_dist <- runif(nrow(MCdata), param[1], param[2])
    } else if (random.regime == "normal") {
      int_dist <- rnorm(nrow(MCdata), param[1], param[2])
    } else if(random.regime == "triangular"){
      int_dist <- EnvStats::rtri(nrow(MCdata), param[1], param[2], param[3])
    } else if (random.regime == "binomial") {
      int_dist <- rbinom(nrow(MCdata), 1, param[1])
    }
    if (incremental == TRUE){
      new_fi <- new_fi + int_dist
    } else {
      new_fi <- int_dist
    }
  }
  new_MCdata <- cbind(MCdata, new_fi)


  return(new_MCdata)
}

ri_pred_surv <- function(BM, MCdatatmp, IT, remove_intercept){

  if(ncol(MCdatatmp)==1) {
    MCdata <- as.data.frame(cbind(1, MCdatatmp))
  } else {
    MCdata <- as.data.frame(MCdatatmp)
  }

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula(
    paste("~", paste(present_vars, collapse = " + "))
  ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]


  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(opts$num_save == 1){
    beta <- BM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(BM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]

  r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) + BM$offset
  eta <- as.numeric( Z %*% beta ) #* BM$sd_Y
  s_hat <- pnorm(r + eta)

  #s_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
  return(s_hat)
}

ri_pred_w <- function(WM, MCdatatmp, IT, remove_intercept) {

  if(ncol(MCdatatmp)==1) {
    MCdata <- as.data.frame(cbind(1, MCdatatmp))
  } else {
    MCdata <- as.data.frame(MCdatatmp)
  }

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula(
    paste("~", paste(present_vars, collapse = " + "))
  ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]

  for(i in 1:ncol(X)) {
    X[,i] <- WM$ecdfs[[i]](X[,i])
  }


  if(opts$num_save == 1){
    beta <- WM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(WM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]


  r <- as.numeric(WM$forest$predict_iteration(as.matrix(X), IT)) + WM$offset
  eta <- as.numeric( Z %*% beta ) #* WM$sd_Y
  w_hat <- pnorm(r + eta)

  #w_hat <- pnorm(as.numeric(WM$forest$predict_iteration(X, IT)) + WM$offset)
  return(w_hat)
}

ri_pred_y <- function(cont, BM, MCdatatmp, IT, remove_intercept) {

  MCdata  <- as.data.frame(MCdatatmp)

  form <- BM$formula
  ln.form <- BM$linear_formula
  opts <- BM$opts

  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  present_vars <- intersect(ln.vars, colnames(MCdata)) # Find which variables are present in the data

  safe_formula <- as.formula( paste("~", paste(present_vars, collapse = " + ") ) ) # Create a new formula with only present variables

  Z <- model.matrix(safe_formula, MCdata)

  if(remove_intercept) {
    idx <- which(colnames(Z) != "(Intercept)")
    Z <- Z[,idx]
  }

  X <- MCdata[, !(names(MCdata) %in% colnames(Z) ) ]


  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(opts$num_save == 1){
    beta <- BM$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(BM$beta[1:opts$num_save,])#/BM$sd_Y
  }

  names(beta) <- ln.vars #### give the names of beta so that the missing beta_k can be found easily at following steps.
  beta <- beta[present_vars]

  if(cont == FALSE) {

    r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) + BM$offset
    eta <- as.numeric(Z %*% beta) #* BM$sd_Y
    y_hat <- pnorm(r + eta)

  } else if(cont == TRUE) {

    r <- as.numeric(BM$forest$predict_iteration(as.matrix(X), IT)) * BM$sd_Y + BM$mu_Y
    eta <- as.numeric(Z %*% beta) #* BM$sd_Y
    y_hat <- r + eta

  }
  return(y_hat)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

##### Check the numebers of prior models and return indices of prior models
### column j, j is larger than the length of variables "X0"(s).
### var.type: a vector introducing the types of variables, marked with "X0", "X", "Fi", "Y", "S", "D" etc..
### time.type: a vector containing the corresponding time point of each variable.
select_history_indices <- function(j, var.type, time.type) {
  vj   <- var.type[j]
  tnum <- as.integer(sub("^T", "", as.character(time.type)))
  tj <- tnum[j] ### current time
  hist_idx <- integer(0)

  if (vj == "S" | vj == "D") {
    for (t in seq_len(tj - 1) ) {
      idx_t <- which(tnum == t & var.type %in% c("S", "D"))
      hist_idx <- c(hist_idx, idx_t)
    }
  }else if (vj == "X") {
    # k = ordinal of this X among X-like's at current time (includes X0 and X)
    idx_now_X <- which(tnum == tj & var.type %in% c("X","X0") )
    k <- match(j, idx_now_X)
    if(is.na(k)) {stop(paste("There is no available X at time", tj)) }

    # In previous time, find the k-th X-like
    for (t in seq_len(tj - 1) ) {
      idx_prev_Xlike <- which(tnum == t & var.type %in% c("X","X0"))
      hist_idx <- c(hist_idx, idx_prev_Xlike[k])
    }
  }else { ## "Y"
    for (t in seq_len(tj - 1) ) {
      idx_t_Y <- which( tnum == t & var.type == "Y")
      hist_idx <- c(hist_idx, idx_t_Y)
    }
  }
  return(hist_idx)
}

##### A function to predict the groups effects when there no historical model
##### if there no missing effect, the function will return the current model
### BM: current models
### linear.formula: the baseline linear formula
### IT: IT-th iteration
pred_beta_BM_noHis <- function(BM, linear_formula, IT){
  opts <- BM$opts

  ln.form <- BM$linear_formula
  ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

  base.ln.vars <- all.vars(linear_formula) # Extract variable names from linear formula

  #Extract terms
  terms_current <- attr(terms(ln.form), "term.labels")
  terms_baseline <- attr(terms(linear_formula), "term.labels")

  # Compare sorted terms
  if (setequal(terms_current, terms_baseline)) {
    BM$linear_formula <- BM$linear_formula
    BM$beta <- BM$beta
  } else {

    beta <- BM$beta#[1:opts$num_save,]
    colnames(beta) <- ln.vars

    not_in_formula <- setdiff(base.ln.vars, ln.vars) # Find names of beta(s) at baseline that are NOT in the current beta
    common_vars <- intersect(ln.vars, base.ln.vars)

    new.beta <- matrix(NA,ncol = length(base.ln.vars), nrow = 1 )
    colnames(new.beta) <- base.ln.vars
    new.beta[,common_vars] <- beta
    BM$beta <- new.beta

    for (rth in 1:opts$num_save) {

      sig_beta <- sd(new.beta, na.rm = TRUE) ### estimate the standard deviation of groups effects
      mean_beta <- mean(new.beta, na.rm = TRUE) ### estimate the mean of groups effects, for genearl bart model, the mean not equals to 0
      pred.beta <- rnorm(length(base.ln.vars)-length(ln.vars), mean = mean_beta, sig_beta)

      BM$beta[rth,not_in_formula] <- pred.beta
    }
    BM$linear_formula <- linear_formula
  }
  return(BM)
}


##### A function to predict the groups effects when there is(are) 1 or more historical models
##### if there no missing effect, the function will return the current model
### BM: current models
### his_BM: a list containing the corresponding historical model(s)
### linear.formula: the baseline linear formula
### IT: IT-th iteration
pred_beta_BM <- function(BM, his_BM, linear_formula, IT, impute_method = NULL, ridge = ridge){

  length_his_BM <- length(his_BM)

   if(length_his_BM == 1){

    opts <- BM$opts
    ln.form <- BM$linear_formula
    ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

    base.ln.vars <- all.vars(linear_formula) # Extract variable names from linear formula at baseline

    #Extract terms
    current_terms <- attr(terms(ln.form), "term.labels")
    base_terms <- attr(terms(linear_formula), "term.labels")

  # Compare sorted terms
    if (setequal(current_terms, base_terms)) {
      BM$linear_formula <- BM$linear_formula
      BM$beta <- BM$beta
    } else {

      his_BM <- his_BM[[1]]
      new.beta <- matrix(NA, nrow = opts$num_save, ncol = length(base.ln.vars) )
      colnames(new.beta) <- base.ln.vars

      for (rth in 1:opts$num_save) {

        #beta <- BM$beta[rth,]
        beta <- as.matrix(BM$beta)[rth, , drop = FALSE]
        his.beta <- his_BM$beta[rth,, drop = FALSE]

        colnames(beta) <- ln.vars
        colnames(his.beta) <- base.ln.vars

        not_in_formula <- setdiff(base.ln.vars, ln.vars) # Find names of beta(s) at baseline that are NOT in the current beta
        common_vars <- intersect(ln.vars, base.ln.vars)

        new.beta[rth, common_vars] <- beta

        full_data.beta  <- t(rbind(his.beta, new.beta[rth,])) %>% as.data.frame()
        colnames(full_data.beta) <- c(paste0("his.beta.", 1:length_his_BM),"beta")
        cc_data.beta <- full_data.beta[complete.cases(full_data.beta),]

        # predict missed beta
        if(impute_method == "BART"){
          mod.beta <- gc_softbart_regression(X = cbind(1, cc_data.beta$his.beta),
                                             Y= cc_data.beta$beta, opts = opts)
          new.X <-  full_data.beta[,!colnames(full_data.beta) %in% "beta"]
          pred.beta <- as.numeric(mod.beta$forest$predict_iteration(as.matrix(cbind(1, new.X)), IT)) * mod.beta$sd_Y +  mod.beta$mu_Y #### replacing all betas
          new.beta[rth, not_in_formula] <- pred.beta[ rownames(full_data.beta) %in% not_in_formula]
        } else if(impute_method == "MVN"){

            full_data.beta <- pred_beta_MVN(
                       full_data.beta = full_data.beta,
                       random = FALSE,
                       ridge = ridge)

          new.beta[rth, ] <- full_data.beta$beta[match(base.ln.vars, rownames(full_data.beta))]
        }

      }

      BM$beta <- new.beta
      BM$linear_formula <- linear_formula
    }

  } else { ### there are more than 1 historical models

    opts <- BM$opts
    ln.form <- BM$linear_formula
    ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

    base.ln.vars <- all.vars(linear_formula) # Extract variable names from linear formula at baseline

    #Extract terms
    current_terms <- attr(terms(ln.form), "term.labels")
    base_terms <- attr(terms(linear_formula), "term.labels")

    # Compare sorted terms
    if (setequal(current_terms, base_terms)) {
      BM$linear_formula <- BM$linear_formula
      BM$beta <- BM$beta
    } else {

      new.beta <- matrix(NA, nrow = opts$num_save, ncol = length(base.ln.vars) )
      colnames(new.beta) <- base.ln.vars

      for (rth in 1:opts$num_save) {

        beta <- as.matrix(BM$beta)[rth, , drop = FALSE]
        colnames(beta) <- ln.vars

        his.beta <- c()
        for (n_prev_idx in 1:length_his_BM) {
          hist_BM_npi <- his_BM[[n_prev_idx]]$beta[rth,]
          his.beta <- rbind(his.beta, hist_BM_npi)
        }

        colnames(his.beta) <- base.ln.vars

        not_in_formula <- setdiff(base.ln.vars, ln.vars) # Find names of beta(s) at baseline that are NOT in the current beta
        common_vars <- intersect(ln.vars, base.ln.vars)

        new.beta[rth,common_vars] <- beta

        full_data.beta  <- t(rbind(his.beta, new.beta[rth,])) %>% as.data.frame()
        colnames(full_data.beta) <- c(paste0("his.beta.", 1:length_his_BM),"beta")
        cc_data.beta <- full_data.beta[complete.cases(full_data.beta),]

        # predict missed beta

        if(impute_method == "BART"){

          mod.beta <- gc_softbart_regression(X = as.matrix(cc_data.beta[,!colnames(cc_data.beta) %in% "beta"]),
                                             Y= cc_data.beta$beta, opts = opts)

          new.X <-  full_data.beta[,!colnames(full_data.beta) %in% "beta"]
          pred.beta <- as.numeric(mod.beta$forest$predict_iteration(as.matrix(new.X), IT)) * mod.beta$sd_Y +  mod.beta$mu_Y #### replacing all betas
          new.beta[rth,not_in_formula] <- pred.beta[ rownames(full_data.beta) %in% not_in_formula]

        } else if(impute_method == "MVN"){

          full_data.beta <- pred_beta_MVN(
            full_data.beta = full_data.beta,
            random = FALSE,
            ridge = ridge)

          new.beta[rth, ] <- full_data.beta$beta[match(base.ln.vars, rownames(full_data.beta))]
        }
      }

      BM$beta <- new.beta
      BM$linear_formula <- linear_formula
    }
  }
  return(BM)
}


#### A function of MVN-imputation
pred_beta_MVN <- function(full_data.beta, random = TRUE, ridge = 1e-6) {

  cc_data.beta <- full_data.beta[complete.cases(full_data.beta), , drop = FALSE]

  mu_hat <- colMeans(cc_data.beta)
  Sigma_hat <- cov(cc_data.beta) + diag(ridge, ncol(full_data.beta))

  obs_idx <- 1:(ncol(full_data.beta) - 1)
  mis_idx <- ncol(full_data.beta)

  mu_o <- mu_hat[obs_idx]
  mu_m <- mu_hat[mis_idx]

  Sigma_oo <- Sigma_hat[obs_idx, obs_idx, drop = FALSE]
  Sigma_mo <- Sigma_hat[mis_idx, obs_idx, drop = FALSE]
  Sigma_mm <- Sigma_hat[mis_idx, mis_idx, drop = FALSE]

  cond_var <- Sigma_mm - Sigma_mo %*% solve(Sigma_oo) %*% t(Sigma_mo)
  cond_var <- max(as.numeric(cond_var), ridge)
  cond_sd <- sqrt(cond_var)

  miss_idx <- which(is.na(full_data.beta[, mis_idx]))

  for (row_i in miss_idx) {

    beta_obs_i <- as.numeric(full_data.beta[row_i, obs_idx])

    cond_mean_i <- mu_m +
      Sigma_mo %*% solve(Sigma_oo) %*% (beta_obs_i - mu_o)

    if (random) {
      full_data.beta[row_i, mis_idx] <- rnorm(1, as.numeric(cond_mean_i), cond_sd)
    } else {
      full_data.beta[row_i, mis_idx] <- as.numeric(cond_mean_i)
    }
  }

  return(full_data.beta)
}

## A function to select subpopulation under strict/weak overlap asumption
## which was apply to the pgs-plugin-gcomp method

select_by_overlap <- function(
  dat,
  expo_cols,
  alive_col,
  group_cols,
  given.regime,
  ovl_typ     = c("weak", "strict"),
  pgs_col     = NULL, #"hat_pgs",
  alive_value = 1
) {
  ovl_typ <- match.arg(ovl_typ, choices = c("weak", "strict"))

  # 1) label each row with a regime id based on exposures
  expo_key <- apply(dat[, expo_cols, drop = FALSE], 1, paste, collapse = "_")
  regime_keys <- vapply(given.regime, paste, collapse = "_", FUN.VALUE = character(1))
  names(regime_keys) <- paste0("regime_", seq_along(regime_keys))

  dat$.regime_id <- NA_character_
  for (i in seq_along(regime_keys)) {
    dat$.regime_id[expo_key == regime_keys[i]] <- names(regime_keys)[i]
  }

  # 2) targeted: in any given regime AND alive
  targeted <- dat[!is.na(dat$.regime_id) & dat[[alive_col]] == alive_value, , drop = FALSE]
  if (nrow(targeted) == 0) return(targeted)

  # helper: for one group, compute overlap across regimes
  group_weak_overlap <- function(g) {
    mins <- numeric(0)
    maxs <- numeric(0)

    for (rid in names(regime_keys)) {
      df_gr <- targeted[targeted[[g]] == 1 & targeted$.regime_id == rid, , drop = FALSE]
      x <- df_gr[[pgs_col]]
      x <- x[!is.na(x)]
      if (length(x) == 0) next

      mins <- c(mins, min(x))
      maxs <- c(maxs, max(x))
    }

    c(weak_min = max(mins), weak_max = min(maxs))
  }

  # 3) weak overlap per group (across regimes)
  weak_by_group <- lapply(group_cols, group_weak_overlap)
  names(weak_by_group) <- group_cols
  weak_min_vec <- vapply(weak_by_group, `[[`, 0.0, "weak_min")
  weak_max_vec <- vapply(weak_by_group, `[[`, 0.0, "weak_max")

  # 4) strict overlap across groups (intersection of group-weak overlaps)
  strict_range <- c(
    strict_min = max(weak_min_vec),
    strict_max = min(weak_max_vec)
  )

  # 5) selection
  in_range <- function(x, lo, hi) !is.na(x) & !is.na(lo) & !is.na(hi) & x >= lo & x <= hi

  if (ovl_typ == "strict") {

    selected <- targeted[
      in_range(targeted[[pgs_col]],
               strict_range["strict_min"],
               strict_range["strict_max"]),
      , drop = FALSE
    ]

  } else {  # ovl_typ == "weak" (each person belongs to exactly one group)

    keep <- rep(FALSE, nrow(targeted))

    for (g in group_cols) {
      idx <- targeted[[g]] == 1
      keep[idx] <- in_range(
        targeted[[pgs_col]],
        weak_by_group[[g]]["weak_min"],
        weak_by_group[[g]]["weak_max"]
      )
    }

    selected <- targeted[keep, , drop = FALSE]
  }

  attr(selected, "strict_overlap_all_groups") <- strict_range
  attr(selected, "weak_overlap_by_group") <- weak_by_group
  selected
}
