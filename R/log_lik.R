#' Compute marginal information criteria for some bayesian linear models
#'
#' @param model The model result (stanfit, stanreg)
#' @param data the input data used to fit the model (if not included in the model
#' object.
#' @param marginalize the random effect to marginalize over
#' @param quadrature_nodes the number of nodes to integrate over
#' @param ... extra parameters for methods
#' @export
compute_marg_lik <- function(model, data, marginalize, quadrature_nodes = 11, ...)
  UseMethod("compute_marg_lik")


#' @export
compute_marg_lik.np_fit <-
  function(model, ...){
    data <- model$fitted[[1]]$data
    ll <- 
      lapply(model$fitted, log_lik) %>%
      Reduce(cbind, .) %>%
      over(dim_l, d ~ c(d[1], d[2] %/% nrow(data), nrow(data))) %>%
      apply(c(1,3), logSumExp)    
    t(ll)
  }

#' @export
compute_marg_lik.edt_fit <-
  function(model
         , marginalize = NULL
         , quadrature_nodes = 11
         , ...
           ){

    . <- NULL
    mutate <- dplyr::mutate
    unnest <- tidyr::unnest
    if(is.null(marginalize)) marginalize <- "ranint_matrix"
      
    data <- model$data
    model <- model$model
    clusters <- as.numeric(data[[marginalize]])

    y <- data$y
    fit <- as.matrix(model, "y_pred")
    ranints <- as.matrix(model, "ranints")
    hsd <- as.matrix(model, "sigma_ranints")
    sigma <- as.matrix(model, "sigma_model")    
    
    if(max(clusters) > ncol(ranints))       
      stop("Some clusters in the model data can't be found in the betas",
           " or appear more than once. This could mean that you tried to",
           " marginalize over something that isn't a random intercept term")
        
    partial_linpred <-
      t(fit - ranints[,clusters])
  
    cluster_frame <-
      data.frame(cluster = unique(clusters)
               , mean = colMeans(ranints)
               , sd = apply(ranints, 2, sd)
               , stringsAsFactors = FALSE)

    std_quad <- gauss.quad.prob(quadrature_nodes, "normal", mu = 0, sigma = 1)
    std_log_weights <- log(std_quad$weights)

    cluster_frame_nodes <-
      cluster_frame %>%
      mutate(node_location = mapply(function(m, s) m + s * std_quad$nodes, mean, sd
                                  , SIMPLIFY = FALSE)
           , pre_weight =
               lapply(sd, function(s, l)
                 log(sqrt(2*pi)) + log(s) + std_quad$nodes^2/2 + std_log_weights)
             ) %>%
      unnest %>%
      mutate(log_weight =
               mapply(function(l,w) w + dnorm(l, sd = hsd, log = TRUE)
                    , node_location, pre_weight
                    , SIMPLIFY = FALSE)) %>%
      mutate(log_lik =
               mapply(function(c,l,w){
                 clusti <- clusters == c
                 yc <- y[clusti]
                 lin_pred <-
                   partial_linpred[clusti,] +
                   l                   

                 ll <-
                   vapply(seq_len(ncol(lin_pred)), function(i){
                     dnorm(yc - lin_pred[,i], sd = sigma[i], log = TRUE)
                   }, FUN.VALUE = numeric(nrow(lin_pred)))

                 t(ll) + as.numeric(w)
               }, cluster, node_location, log_weight,
               SIMPLIFY = FALSE)) %>%
       group_by(cluster) %>%
       summarize(log_lik =
         list(
           Reduce(function(acc, x){                  
             acc <- c(acc,x)
             dim(acc) <- c(dim(x), length(acc) %/% prod(dim(x)))
             acc
           }, log_lik) %>%
           apply(1, logSumExp) #sum along 2 (structures) and 3 (nodes)
         ))

    Reduce(rbind, cluster_frame_nodes$log_lik) %>%
      `rownames<-`(cluster_frame_nodes$cluster)
  }

#' @export
compute_marg_lik.stanreg <-
  function(model
         ,  marginalize = NULL
         ,  quadrature_nodes = 11
         , ...
           ){

    . <- NULL
    mutate <- dplyr::mutate
    unnest <- tidyr::unnest

    data <- model$data
    if(is.null(marginalize)) marginalize <- "ID"
    clusters <- as.character(unique(data[[marginalize]]))
    args <- rstanarm:::ll_args.stanreg(model)
    xdata <- as.matrix(rstanarm:::.xdata(args$data))
    y <- args$data$y
    
    if(!all(clusters %in% colnames(xdata)) ||
       ncol(xdata[,clusters]) != length(clusters))       
      stop("Some clusters in the model data can't be found in the betas",
           " or appear more than once. This could mean that you tried to",
           " marginalize over something that isn't a random intercept term")

    beta <- args$draws$beta
    which_bs <- match(clusters, colnames(xdata))
    b_samples <- args$draws$beta[,which_bs]
    hvar <- as.data.frame(model)[, paste0("Sigma[", marginalize, ":(Intercept),(Intercept)]")]       
        
    partial_linpred <-
      tcrossprod(xdata[,-which_bs, drop = FALSE]
               , beta[,-which_bs, drop = FALSE]) %>%
      ## gpuR::tcrossprod(
      ##         vclMatrix(xdata[,-which_bs, drop = FALSE])
      ##       , vclMatrix(beta[,-which_bs, drop = FALSE])) %>%
      as.matrix
          
    cluster_frame <-
      data.frame(cluster = clusters
               , b = which_bs
               , bnm = colnames(beta)[which_bs]
               , mean = colMeans(b_samples)
               , sd = apply(b_samples, 2, sd)
               , stringsAsFactors = FALSE)

    std_quad <- gauss.quad.prob(quadrature_nodes, "normal", mu = 0, sigma = 1)
    std_log_weights <- log(std_quad$weights)

    cluster_frame_nodes <-
      cluster_frame %>%
      mutate(node_location = mapply(function(m, s) m + s * std_quad$nodes, mean, sd
                                  , SIMPLIFY = FALSE)
           , pre_weight =
               lapply(sd, function(s, l)
                 log(sqrt(2*pi)) + log(s) + std_quad$nodes^2/2 + std_log_weights)
             ) %>%
      unnest %>%
      mutate(log_weight =
               mapply(function(l,w) w + dnorm(l, sd = sqrt(hvar), log = TRUE)
                    , node_location, pre_weight
                    , SIMPLIFY = FALSE)) %>%
      mutate(log_lik =
               mapply(function(c,l,w,b){
                 clusti <- args$data[,c] == 1
                 yc <- y[clusti]
                 lin_pred <-
                   partial_linpred[clusti,] +
                   xdata[clusti, b] * l                   

                 ll <-
                   vapply(seq_len(ncol(lin_pred)), function(i){
                     dnorm(yc - lin_pred[,i], sd = args$draws$sigma[i], log = TRUE)
                   }, FUN.VALUE = numeric(nrow(lin_pred)))

                 t(ll) + w
               }, cluster, node_location, log_weight, b,
               SIMPLIFY = FALSE)) %>%
       group_by(cluster) %>%
       summarize(log_lik =
         list(
           Reduce(function(acc, x){                  
             acc <- c(acc,x)
             dim(acc) <- c(dim(x), length(acc) %/% prod(dim(x)))
             acc
           }, log_lik) %>%
           apply(1, logSumExp) #sum along 2 (structures) and 3 (nodes)
         ))

    Reduce(rbind, cluster_frame_nodes$log_lik) %>%
      `rownames<-`(cluster_frame_nodes$cluster)
  }

dont <- function(...){invisible(NULL)}

dont({
## Test marginal_likelihood
d <-
  data_frame(ID = 1:10, group = sample(0:1, 10, replace = TRUE)) %>%
  mutate(y = map(ID, ~ (ID - 5)/10 + (group - .5)*rnorm(10) + rnorm(10))) %>%
  unnest

m <- stan_lmer(y ~ group + (1 | ID), data = d)

mar <- compute_log_lik(m, "ID")

hs <- as.matrix(m$stanfit)[,"Sigma[ID:(Intercept),(Intercept)]"]
id1s <- as.matrix(m$stanfit)[,"b[(Intercept) ID:1]"]
s <- as.matrix(m$stanfit)[,"sigma"]

std_quad <- gauss.quad.prob(11, "normal", mu = 0, sigma = 1)
std_log_weights <- log(std_quad$weights)
sites <- mean(id1s) + sd(id1s) * std_quad$nodes
weights <-
  vapply(hs, function(hsi){
    log(sqrt(2*pi)) + log(sd(id1s)) + std_quad$nodes^2/2 +
      dnorm(sites, sd = sqrt(hsi), log = TRUE) + std_log_weights
  }, numeric(length(sites)))

args <- rstanarm:::ll_args.stanreg(m)
linp <- tcrossprod(args$draws$beta, as.matrix(rstanarm:::.xdata(args$data)))
linp <- linp - id1s ## remove indiv effect
linp <- linp[,d$ID == 1]
yi <- d$y[d$ID == 1]

ll <- sapply(seq_len(length(id1s)), function(i){
    sapply(seq_len(nrow(weights)), function(j){  
      dnorm(linp[i,] + sites[j] - yi, sd = s[i], log = TRUE) + weights[j,i]
    })
}) %>%
  apply(2, logSumExp)

## Do they match?
all.equal(mar[1,] , ll)
})
