#' Get the fitted values from bayesian hierarchical linear models
#'
#' Much faster for gaussian linear models than the rstanarm equivalent
#'
#' @param mod The model
#' @export
get_fitted <- function(mod) UseMethod("get_fitted")

#' @export
get_fitted.edt_fit <- function(mod) t(as.matrix(mod$model, "y_pred"))

#' @export
get_fitted.stanreg <- function(mod){
    args <- rstanarm:::ll_args.stanreg(mod)
    xdata <- as.matrix(rstanarm:::.xdata(args$data))
    beta <- args$draws$beta
    tcrossprod(xdata, beta)
}

#' @export
get_fitted.np_fit <- function(mod){
  lapply(mod$fitted, get_fitted) %>%
    reduce(rbind)
}

#' @export
get_fitted.npl_fit <- function(mod){
  lapply(mod$fitted, get_fitted) %>%
    reduce(rbind)
}

#' Get posterior predictions from a bayesian hierarchical linear model
#'
#' Faster than the rstanarm version.
#'
#' @param mod The model
#' @export
post_pred <- function(mod) UseMethod("post_pred")

#' @export
post_pred.edt_fit <- function(mod){
  fit <- get_fitted(mod)
  err <-
    as.matrix(mod$model, "sigma_model")[,1] %>%
    lapply(function(s) rnorm(nrow(fit), sd = s)) %>%
    reduce(cbind)

  fit + err
}

#' @export
post_pred.stanreg <- function(mod){
  fit <- get_fitted(mod)
  err <- as.matrix(mod$stanfit, "aux")[,1] %>%
    lapply(function(s) rnorm(nrow(fit), sd = s)) %>%
    reduce(cbind)

  fit + err
}

#' @export
post_pred.np_fit <- function(mod){
  lapply(mod$fitted, post_pred) %>%
    reduce(rbind)
}

#' @export
post_pred.npl_fit <- function(mod){
  lapply(mod$fitted, post_pred) %>%
    reduce(rbind)
}

#' Get the y value for a model, useful in the case of no-pooling
#' where the variables get shuffled.
#' @param mod The model of interest
#' @export
extract_y <- function(mod) UseMethod("extract_y")

#' @export
extract_y.edt_fit <- function(mod){
  mod$data$y
}

#' @export
extract_y.np_fit <- function(mod){
  sapply(mod$fitted, extract_y)        
}

#' @export
extract_y.npl_fit <- function(mod){
  sapply(mod$fitted, extract_y)        
}

#' @export
extract_y.stanreg <- function(mod){
  mod$y  
}

#' Get the data for a model, useful in the case of no-pooling
#' where the variables get shuffled.
#' @param mod The model of interest
#' @export
extract_data <- function(mod) UseMethod("extract_data")

#' @export
extract_data.edt_fit <- function(mod) mod$data

#' @export
extract_data.stanreg <- function(mod) mod$data

#' @export
extract_data.np_fit <- function(mod) bind_rows(lapply(mod$fitted, extract_data))

#' @export
extract_data.npl_fit <- function(mod) bind_rows(lapply(mod$fitted, extract_data))


#' Compute subject specific sum of squares
#'
#' @param The predicted value
#' @param The model of interest
#' @export
compute_subject_error <- function(pred, mod){
  y <- extract_y(mod)
  mod_d <- extract_data(mod)

  if(inherits(mod, "edt_fit")){
    id <- mod_d$ranint_matrix[,1]
  } else {
    id <- mod_d$ID
  }
  
  tapply(seq_len(nrow(pred)), list(id), function(rows){
    colSums((y[rows] - pred[rows,,drop = FALSE])^2)
  }) %>%
  reduce(rbind) %>%
  t
}
