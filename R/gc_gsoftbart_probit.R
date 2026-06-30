gc_gsoftbart_probit <- function(X, Y, Z, num_tree = 20,
                            k = 1, hypers = NULL, opts = NULL, verbose = TRUE) {

  ## Get design matricies and groups for categorical

  # dv <- dummyVars(formula, data)
  # terms <- attr(dv$terms, "term.labels")
  # group <- dummy_assign(dv)
  # suppressWarnings({
  #   X_train <- predict(dv, data)
  #   X_test  <- predict(dv, test_data)
  # })
  # Y_train <- model.response(model.frame(formula, data))
  # Y_test  <- model.response(model.frame(formula, test_data))

  Y_train <- Y #model.response(model.frame(formula, data))
  X_train <- as.matrix(X) ## covariates in non-linear function
  Z_train <- as.matrix(Z) ## covariates in linear function

  if( length(colnames(Z)) != 0 ){
      linear_formula <- as.formula(paste("~", paste0(colnames(Z), collapse  = "+")))
  } else { linear_formula <- as.formula(paste("~", paste0("Z.", 1:ncol(Z),collapse  = "+"))) }

  stopifnot(length(table(Y_train)) == 2)

  pnorm_offset <- mean(Y_train)
  offset <- qnorm(pnorm_offset)

  ## Set up hypers

  if(is.null(hypers)) {
    hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train, tgroup = (1:ncol(X_train) - 1) )
  } else {
    hypers$X <- cbind(X_train, Z_train)
    hypers$Y <- Y_train
  }

  hypers$sigma_mu = 3 / k / sqrt(num_tree)
  hypers$sigma <- 1
  hypers$sigma_hat <- 1
  hypers$num_tree <- num_tree
  hypers$group <- (1:ncol(X_train) - 1)

  ## Set up opts

  if(is.null(opts)) {
    opts <- Opts()
  }
  opts$update_sigma <- FALSE
  opts$num_print <- 2147483647

  ## Normalize!

  make_01_norm <- function(x) {
    a <- min(x)
    b <- max(x)
    return(function(y) (y - a) / (b - a))
  }

  ecdfs   <- list()
  for(i in 1:ncol(X_train)) {
    ecdfs[[i]] <- ecdf(X_train[,i])
    if(length(unique(X_train[,i])) == 1) ecdfs[[i]] <- identity
    if(length(unique(X_train[,i])) == 2) ecdfs[[i]] <- make_01_norm(X_train[,i])
  }

  for(i in 1:ncol(X_train)) {
    X_train[,i] <- ecdfs[[i]](X_train[,i])
  }

  ## Make forest ----

  probit_forest <- MakeForest(hypers, opts)

  ## Initialize output

  r_train      <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  eta_train    <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  mu_train     <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  beta_out     <- matrix(NA, nrow = opts$num_save, ncol = ncol(Z_train))
  #sigma_out    <- numeric(opts$num_save)
  sigma_mu_out <- numeric(opts$num_save)
  varcounts <- matrix(NA, nrow = opts$num_save, ncol = ncol(X_train))

  ## Prepare running the chain ----
  r     <- as.numeric(probit_forest$do_predict(X_train))
  beta  <- numeric(ncol(Z_train))
  sigma <- 1 #forest$get_sigma() ## here sigma always equal to 1
  eta   <- numeric(nrow(Z_train))

  mu <- r + eta + offset
  lower <- ifelse(Y_train == 0, -Inf, 0)
  upper <- ifelse(Y_train == 0, 0, Inf)

  ## Function for updating beta ----
  ZtZi <- solve(t(Z_train) %*% Z_train)

  update_beta <- function(R, Z, V, sigma) {
    beta_hat <- as.numeric(V %*% t(Z) %*% R)
    beta <- MASS::mvrnorm(n = 1, mu = beta_hat, Sigma = sigma^2 * V)
    return(as.numeric(beta))
  }

  ## Warmup
  for(i in 1:opts$num_burn) {

        W <- truncnorm::rtruncnorm(n = length(Y_train), a = lower, b = upper,
                      mean = mu, sd = 1)

      ## Update beta ----
      R <- W - r - offset
      beta <- update_beta(R, Z_train, ZtZi, sigma^2)
      eta <- as.numeric(Z_train %*% beta)

      ## Update forest ----
      R <- W - eta - offset
      r <- as.numeric(probit_forest$do_gibbs(X_train, R , X_train, 1))

      ## Update mu
      mu <- r + eta + offset
       }

  ## Save
  for(i in 1:opts$num_save) {
    for(j in 1:opts$num_thin) {
        ## Sample W
        W <- truncnorm::rtruncnorm(n = length(Y_train), a = lower, b = upper,
                        mean = mu, sd = 1)

        ## Update beta ----
        R <- W - r - offset
        beta <- update_beta(R, Z_train, ZtZi, sigma^2)
        eta <- as.numeric(Z_train %*% beta)

        ## Update forest ----
        R <- W - eta - offset
        r <- as.numeric(probit_forest$do_gibbs(X_train, R, X_train, 1))

        ## update mu  ----
        mu <- r + eta + offset

    }

    r_train[i,] <- as.numeric(probit_forest$do_predict(X_train))
    eta_train[i,] <- as.numeric(Z_train %*% beta)
    mu_train[i,] <- r_train[i,] + eta_train[i,] + offset
    beta_out[i,] <- beta
    sigma_mu_out[i] <- probit_forest$get_sigma_mu()
    varcounts[i,] <- as.numeric(probit_forest$get_counts())


  }

  p_train <- pnorm(mu_train)

  colnames(varcounts) <- colnames(X_train)

  out <- list(
   r_train = r_train, eta_train = eta_train, mu_train = mu_train,
    var_counts = varcounts,
    beta = beta_out,
    sigma_mu = sigma_mu_out,
    #mu_train_mean = colMeans(mu_train),
    #p_train_mean = colMeans(p_train),
    offset = offset,
    pnorm_offset = pnorm(offset),
    p_train = p_train,
    #formula = formula,
    ecdfs = ecdfs,
    opts = opts,
    forest = probit_forest,
    linear_formula = linear_formula)

  class(out) <- "gc_gsoftbart_probit"
  return(out)

}
