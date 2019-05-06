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
