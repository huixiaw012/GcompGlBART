#' @importFrom stats as.formula complete.cases cov ecdf lm model.matrix na.omit pnorm predict qnorm quantile rbinom rmultinom rnorm runif sd terms
#' @importFrom GcompBART Hypers Opts MakeForest gc_softbart_regression
#' @importFrom magrittr %>%
NULL

utils::globalVariables(c("BM", "WModels", "num_save", "num_thin"))
