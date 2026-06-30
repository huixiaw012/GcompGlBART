#### Z should have colnames
gc_gsoftbart_regression <- function(X, Y, Z, num_tree = 20, k = 2,
                                    hypers = NULL, opts = NULL,
                                    verbose = TRUE, warn = TRUE ) {

  ## Get design matricies and groups for categorical
  # data <- data.frame(X = X, Y = Y)
  # dv <- dummyVars(Y ~ X)
  # terms <- attr(dv$terms, "term.labels")
  # group <- dummy_assign(dv)
  # suppressWarnings({
  #   X_train <- predict(dv, data)
  # })

  Y_train <- Y #model.response(model.frame(formula, data))
  X_train <- as.matrix(X) ## covariates in non-linear function
  Z_train <- as.matrix(Z) ## covariates in linear function


  if( length(colnames(Z)) != 0 ){
      linear_formula <- as.formula(paste("~", paste0(colnames(Z), collapse  = "+")))
  } else { linear_formula <- as.formula(paste("~", paste0("Z.", 1:ncol(Z),collapse  = "+"))) }

  stopifnot(is.numeric(Y_train))

  mu_Y <- mean(Y_train)
  sd_Y <- sd(Y_train)
  Y_train <- (Y_train - mu_Y) / sd_Y

  ## Set up hypers
  if(is.null(hypers)) {
    hypers <- Hypers(X = cbind(X_train, Z_train), Y = Y_train, normalize_Y = FALSE,
                     tgroup = (1:ncol(X_train) - 1)) #, tgroup_size = table(1:ncol(X_train) - 1))
  } else {
    hypers$X <- cbind(X_train, Z_train)
    hypers$Y <- Y_train
  }
  hypers$sigma_mu = 3 / k / sqrt(num_tree)
  hypers$num_tree <- num_tree
  hypers$group <- (1:ncol(X_train) - 1) #group. !!!!!! here do we need to change something?

  ## Set up opts

  if(is.null(opts)) {
    opts <- Opts()
  }
  opts$num_print <- .Machine$integer.max

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
  reg_forest <- MakeForest(hypers, opts)

  ## Initialize output

  r_train      <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  eta_train    <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  mu_train     <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
  beta_out     <- matrix(NA, nrow = opts$num_save, ncol = ncol(Z_train))
  sigma_out    <- numeric(opts$num_save)
  sigma_mu_out <- numeric(opts$num_save)
  varcounts <- matrix(NA, nrow = opts$num_save, ncol = ncol(X_train))

  ## Prepare running the chain ----
  r     <- as.numeric(reg_forest$do_predict(X_train))
  beta  <- numeric(ncol(Z_train))
  sigma <- reg_forest$get_sigma()
  eta   <- numeric(nrow(Z_train))

  ## Function for updating beta ----
  ZtZi <- solve(t(Z_train) %*% Z_train)

  update_beta <- function(R, Z, V, sigma) {
    beta_hat <- as.numeric(V %*% t(Z) %*% R)
    beta <- MASS::mvrnorm(n = 1, mu = beta_hat, Sigma = sigma^2 * V)
    return(as.numeric(beta))
  }


  ## Warmup

  for(i in 1:opts$num_burn) {

    ## Update beta ----
    R <- Y_train - r
    beta <- update_beta(R, Z_train, ZtZi, sigma^2)
    eta <- as.numeric(Z_train %*% beta)

    ## Update forest and sigma ----
    R <- Y_train - eta
    r <- as.numeric(reg_forest$do_gibbs(X_train, R, X_train, 1))
    sigma <- reg_forest$get_sigma()

  }

  ## Save
  for(i in 1:opts$num_save) {
    for(j in 1:opts$num_thin) {

       ## Update beta ----
      R <- Y_train - r
      beta <- update_beta(R, Z_train, ZtZi, sigma^2)
      eta <- as.numeric(Z_train %*% beta)

      ## Update forest and sigma ----
      R <- Y_train - eta
      r <- as.numeric(reg_forest$do_gibbs(X_train, R, X_train, 1))
      sigma <- reg_forest$get_sigma()
    }

    r_train[i,] <- as.numeric(reg_forest$do_predict(X_train)) * sd_Y + mu_Y
    eta_train[i,] <- as.numeric(Z_train %*% beta) * sd_Y
    mu_train[i,] <- r_train[i,] + eta_train[i,]
    beta_out[i,] <- beta * sd_Y
    sigma_out[i] <- sigma * sd_Y
    sigma_mu_out[i] <- reg_forest$get_sigma_mu() * sd_Y
    varcounts[i,] <- as.numeric(reg_forest$get_counts())

  }

  colnames(varcounts) <- colnames(X_train)

  out <- list(
    r_train = r_train, eta_train = eta_train, mu_train = mu_train,
    beta = beta_out, sigma = sigma_out,
    #sigma_mu = sigma_mu_out,
    var_counts = varcounts,
    #mu_train = mu_train,
    #mu_train_mean = colMeans(mu_train),
    opts = opts, ecdfs = ecdfs,
    mu_Y = mu_Y, sd_Y = sd_Y,
    forest = reg_forest,
    linear_formula = linear_formula)

  class(out) <- "gc_gsoftbart_regression"

  return(out)
}
