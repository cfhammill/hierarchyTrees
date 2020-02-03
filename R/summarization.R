#' Summarize a bayesian linear model
#'
#' Extracts the salient results from a tree model. The new interface to
#' the old "get_*_results" functions
#'
#' @param model the model to summarize
#' @param tree the tree used to fit the model
#' @param leaf_scales scaling factors for the leaves, useful
#' if the leaf response variables were scaled before modelling.
#' @param add_fixed_effect whether or not to incorporate the fixed
#' add the main effect to the leaf effects.
#' @param p0 the leaf coefficient name
#' @param p1 the first layer parent name
#' @param p2 the second layer parent name
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
summarize_model <-
    function(model
           , tree
           , leaf_scales = NULL
           , just_leaves = TRUE
           , ...
             ){
        UseMethod("summarize_model")
    }

#' @export
summarize_model.h0_fit <-
    function(model, tree, leaf_scales = NULL, add_fixed_effect = TRUE
           , main = "group", p0 = "p0", ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_sglm_results(model
                     , tree = tree
                     , sds = leaf_scales
                     , addFix = add_fixed_effect
                     , justLeaves = TRUE
                     , main = main, struct = p0)
    }

#' @export
summarize_model.h1_fit <-
    function(model, tree, leaf_scales = NULL, add_fixed_effect = TRUE
           , main = "group", p0 = "p0", p1 = "p1", ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_hsglm_results(model
                     , tree = tree
                     , sds = leaf_scales
                     , addFix = add_fixed_effect
                     , main = main, struct = p0
                     , parent = p1)
    }


#' @export
summarize_model.h2_fit <-
    function(model, tree, leaf_scales = NULL, add_fixed_effect = TRUE
           , main = "group", p0 = "p0", p1 = "p1", p2 = "p2"
           , ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_h2sglm_results(model
                     , tree = tree
                     , sds = leaf_scales
                     , addFix = add_fixed_effect
                     , main = main, struct = p0
                     , parent = p1, grandparent = p2)
    }

#' @export
summarize_model.edt_fit <-
    function(model, tree, leaf_scales = NULL, add_fixed_effect = TRUE, ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_ept_results(model$model
                     , tree = tree
                     , sds = leaf_scales
                     , addFix = add_fixed_effect
                     , justLeaves = TRUE
                       )
    }

#' @export
summarize_model.np_fit <-
    function(model, tree, leaf_scales = NULL, main = "group", ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_np_results(model
                     , tree = tree
                     , sds = leaf_scales
                     , justLeaves = TRUE
                     , main = main
                       )
    }

#' @export
summarize_model.npl_fit <-
    function(model, tree, leaf_scales = NULL, main = "group", ...){
        if(is.null(leaf_scales)){
            leaves <- tree$Get("name", filterFun = isLeaf)
            leaf_scales <-
                rep(1, length(leaves)) %>%
                setNames(leaves)                
        }
          
        get_npl_results(model
                     , tree = tree
                     , sds = leaf_scales
                     , justLeaves = TRUE
                     , main = main
                       )
    }

#' Extract the salient results from an ept mod
#' 
#' @param ept_mod The ept_mod
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param justLeaves Whether or not to just extract effects at leaves
#' @param addFix Whether to add the fixed effects before scaling
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
get_ept_results <-
  function(ept_mod, tree, sds, justLeaves = FALSE, addFix = TRUE){
    nodes <- tree$Get("name")
    lnodes <- tree$Get("name", filterFun = isLeaf)
    
    post <- extract(ept_mod)$b[ , order(nodes), 2]
    colnames(post) <- sort(nodes)

    if(justLeaves){
      post <- post[,lnodes]
      nodes <- lnodes
    }

    b_fix <-
      extract(ept_mod)$b_fix 

    if(addFix)
      post <- post[,nodes] + b_fix[,2]

    scaled_post <- t(t(post[,nodes]) * sds[nodes])
    
    raw_effects <-
      apply(post[,nodes], 2, median)    

    scaled_effects <-
      apply(scaled_post[,nodes], 2, median)

    ranints <-
      extract(ept_mod) %>%
      .$ranints %>%
      apply(2,median) %>%
      `*`(mean(sds[nodes]))

    list(fix_post = b_fix
       , fix = b_fix %>% apply(2, median)
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param justLeaves Whether or not to just extract effects at leaves
#' @param addFix Whether to add the fixed effects before scaling 
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
get_sglm_results <-
  function(smod, tree, sds, justLeaves = FALSE, addFix = TRUE, main = "SEX", struct = "name"){
    ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
    nodes <- tree$Get("name", filterFun = ff)
    post <- as.matrix(smod)
    cols <- colnames(post)

    effect_name <- function(x){
      paste0("b\\[", main, ".*", struct, ":", x, "]")
    }
    effect_cols <- sapply(nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(effect_name(.), cols) }
    })

    b_post <- post[,effect_cols]
    fix <- post[ , c("(Intercept)",main)] 
    colnames(b_post) <- nodes

    if(addFix)
      b_post <- b_post + fix[,2]

    scaled_post <- t(t(b_post) * (sds))
      
    raw_effects <- 
      b_post %>%
      apply(2, median) 

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- as.numeric(ranef(smod)$ID[,1])
   
    list(fix_post = fix
       , fix = fix %>% apply(2, median)
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param justLeaves Whether or not to just extract effects at leaves
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
get_fp_results <-
  function(smod, tree, sds, justLeaves = FALSE, main = "SEX"){
    ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
    nodes <- tree$Get("name", filterFun = ff)
    post <- as.matrix(smod)
    cols <- colnames(post)
    
    b_post <- post[ , rep(main, length(nodes))]
    fix <- post[ , c("(Intercept)", main)] %>% apply(2, median)

    scaled_post <- t(t(b_post) * (sds))

    raw_effects <- 0

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- as.numeric(ranef(smod)$ID[,1])
   
    list(fix = fix
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param justLeaves Whether or not to just extract effects at leaves
#' @return A list containing
#' 1. `fix` 0 for this model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` 0 for this model
#' 4. `b_post` The posterior for the effects
#' @export
get_np_results <-
  function(smod, tree, sds, justLeaves = FALSE, main = "SEXM"){
    ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
    nodes <- tree$Get("name", filterFun = ff)

    if(!"name" %in% names(smod))
        smod$name <- smod$p0
    
    post <- sapply(nodes, function(n){
      mod <- filter(smod, name == n)$fitted[[1]]
      as.matrix(mod)[ , main, drop = FALSE] %>%
        `colnames<-`(n)
    })
    
    b_post <- post
    
    fix <- 0

    scaled_post <- t(t(b_post) * (sds))

    raw_effects <- 
      b_post %>%
      apply(2, median) 

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- 0
   
    list(fix = fix
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param justLeaves Whether or not to just extract effects at leaves
#' @return A list containing
#' 1. `fix` 0 for this model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` 0 for this model
#' 4. `b_post` The posterior for the effects
#' @export
get_npl_results <-
  function(smod, tree, sds, justLeaves = FALSE, main = "SEXM"){
    ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
    nodes <- tree$Get("name", filterFun = ff)

    if(!"name" %in% names(smod))
        smod$name <- smod$p0

    post <- sapply(nodes, function(n){
      mod <- filter(smod, name == n)$fitted[[1]]
      as.matrix(mod)[ , main, drop = FALSE] %>%
        `colnames<-`(n)
    })
    
    b_post <- post
    
    fix <- 0

    scaled_post <- t(t(b_post) * (sds))

    raw_effects <- 
      b_post %>%
      apply(2, median) 

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- 0
   
    list(fix = fix
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param addFix Whether to add the fixed effects before scaling
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
get_hsglm_results <-
  function(smod, tree, sds, addFix = TRUE, main = "SEX"
         , struct = "name", parent = "parent"){
    leaf_nodes <- tree$Get("name", filterFun = isLeaf)
    parent_nodes <- tree$Get("name", filterFun = Negate(isLeaf))
    
    post <- as.matrix(smod)
    cols <- colnames(post)
    effect_name <- function(x){
      paste0("b\\[", main, ".*", struct, ":", x, "]")
    }
    effect_cols <- sapply(leaf_nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(effect_name(.), cols) }
    })

    peff_name <- function(x){
      paste0("b\\[", main, ".*", parent, ":", x, "]")
    }
    peff_cols <- sapply(parent_nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(peff_name(.), cols) } %>%
        { `if`(length(.) == 0, NA, .) }
    })

    tree$Do(function(n){
      n$post <- post[,effect_cols[n$name]]
    }, filterFun = isLeaf)

    tree$Do(function(n){ 
      n$post <- post[,peff_cols[n$name]] 
    }, filterFun = function(n){
      !isLeaf(n) && any(sapply(n$children, isLeaf))
    })

    tree$Do(function(n){
      if(!is.null(n$parent) && !is.null(n$parent$post))
        n$post <- n$post + n$parent$post
    }, traversal = "pre-order", filterFun = isLeaf)

    b_post <-
      tree$Get("post", filterFun = isLeaf, simplify = FALSE) %>%
      simplify2array
    
    fix <- post[ , c("(Intercept)",main)] 
    colnames(b_post) <- leaf_nodes

    if(addFix)
      b_post <- b_post + fix[,2]

    scaled_post <- t(t(b_post) * (sds))

    raw_effects <- 
      b_post %>%
      apply(2, median) 

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- as.numeric(ranef(smod)$ID[,1])
   
    list(fix_post = fix 
       , fix = fix %>% apply(2, median)
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Extract the salient results from an stan_glmer model
#' 
#' @param smod The [stan_glmer] model
#' @param tree The tree of interest
#' @param sds SDs used in scaling leaf volumes
#' @param addFix Whether to add the fixed effects before scaling
#' @return A list containing
#' 1. `fix` The fixed effects of the model
#' 2. `effects` The estimated median effects in the model
#' 3. `ranints` The estimated random intercepts from the model
#' 4. `b_post` The effect posterior samples
#' @export
get_h2sglm_results <-
  function(smod, tree, sds, addFix = TRUE, main = "SEX"
         , struct = "name", parent = "parent", grandparent = "gparent"){
    is_p1 <- function(n) !isLeaf(n) && any(sapply(n$children, isLeaf))
    
    leaf_nodes <- tree$Get("name", filterFun = isLeaf)
    p1_nodes <- tree$Get("name", filterFun = is_p1)
    p2_nodes <- tree$Get("name", filterFun = function(n){
      !isLeaf(n) && any(sapply(n$children, is_p1))
    })
    
    post <- as.matrix(smod)
    cols <- colnames(post)

    effect_name <- function(x){
      paste0("b\\[", main, ".*", struct, ":", x, "]")
    }
    effect_cols <- sapply(leaf_nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(effect_name(.), cols) }
    })

    peff_name <- function(x){
      paste0("b\\[", main, ".*", parent, ":", x, "]")
    }
    peff_cols <- sapply(p1_nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(peff_name(.), cols) } %>%
        { `if`(length(.) == 0, NA, .) }
    })
    
    geff_name <- function(x){
        paste0("\\[", main, ".*", grandparent, ":", x, "]")
    }
    
    geff_cols <- sapply(p2_nodes, function(n){
      n %>%
        gsub("([().])", "\\\\\\1", .) %>%
        gsub(" ", "_", .) %>% { grep(geff_name(.), cols) } %>%
        { `if`(length(.) == 0, NA, .) }
    })

    tree$Do(function(n){
      n$post <- post[,effect_cols[n$name]]
    }, filterFun = isLeaf)

    tree$Do(function(n){ 
      n$post <- post[,peff_cols[n$name]] 
    }, filterFun = is_p1)

    tree$Do(function(n){ 
      n$post <- post[,geff_cols[n$name]] 
    }, filterFun = function(n){
      !isLeaf(n) && any(sapply(n$children, is_p1))
    })

    tree$Do(function(n){
      if(!is.null(n$parent) && !is.null(n$parent$post))
        n$post <- n$post + n$parent$post
    }, traversal = "pre-order", filterFun = isLeaf)

    b_post <-
      tree$Get("post", filterFun = isLeaf, simplify = FALSE) %>%
      simplify2array
    
    fix <- post[ , c("(Intercept)",main)]
    colnames(b_post) <- leaf_nodes

    if(addFix)
      b_post <- b_post + fix[,2]
    
    scaled_post <- t(t(b_post) * (sds))
    
    raw_effects <- 
      b_post %>%
      apply(2, median) 

    scaled_effects <-
      scaled_post %>% apply(2, median)

    ranints <- as.numeric(ranef(smod)$ID[,1])
   
    list(fix_post = fix
       , fix = fix %>% apply(2, median)
       , effects = scaled_effects
       , raw_effects = raw_effects
       , ranints = ranints
       , b_post = b_post
       , scaled_post = scaled_post)
  }

#' Compute the likelihood for a tree model
#'
#' @param mod A stan model with at least `y_pred`
#' and `sigma_model`
#' @param y The y data
#' @return A matrix with the same dimensions as
#' `extract(mod, "y_pred")[[1]]` i.e. samples x obs
#' containing the gaussian likelihood of each
#' y_pred sample.
#' @export
logLik_ept <- function(mod, y){
  y_pred <- extract(mod, "y_pred")[[1]]
  sigma <- extract(mod, "sigma_model")[[1]]

  t(sapply(seq_len(nrow(y_pred)),
         function(i) log(dnorm(y - y_pred[i,], sd = sigma[i]))))
}

#' Compute the likelihood for a tree model
#'
#' @param mod A data.frame with unpooled models
#' @return A log-likelihood matrix
#' @export
logLik_np <- function(mod){
  lapply(mod$fitted, log_lik) %>%
    Reduce(cbind, .) 
}

#' Compute the marginal likelihood of a data point
#'
#' This assumes no knowledge of the individual random intercept from
#' and simulates unobserved subjects.
#'
#' @param The model
#' @param y the observed data
#' @param indiv an integer index to the individual
#' @return A matrix with the same dimensions as
#' `extract(mod, "y_pred")[[1]]` i.e. samples x obs
#' containing the gaussian likelihood of each
#' y_pred sample.
#' @export
logLik_ept_nocluster <- function(mod, y, indiv){
  y_pred <- extract(mod, "y_pred")[[1]]
  sigma <- extract(mod, "sigma_model")[[1]]
  ranints <- as.matrix(real_edt, pars = "ranints")

  per_y_ranints <- ranints[,indiv]
  y_pred_minus_ranints <- y_pred - per_y_ranints

  t(sapply(seq_len(nrow(y_pred_minus_ranints)),
           function(i) log(dnorm(y - y_pred_minus_ranints[i,]
                               , sd = sigma[i]))))
}

#' Compute the marginal likelihood of a data point
#'
#' This assumes no knowledge of the individual random intercept from
#' and simulates unobserved subjects.
#'
#' @param The model
#' @param y the observed data
#' @param indiv a character vector index to the individual
#' @return A matrix with the same dimensions as
#' `extract(mod, "y_pred")[[1]]` i.e. samples x obs
#' containing the gaussian likelihood of each
#' y_pred sample.
#' @export
logLik_sglm_nocluster <- function(mod, y, indiv){
  y_pred <- posterior_linpred(mod)
  sigma <- extract(mod$stanfit, "aux")[[1]]
  ranints <-
    as.matrix(mod) %>%
    { .[ , grepl("b\\[\\(Intercept\\) ID.*", colnames(.))] }

  indiv <- paste0("b[(Intercept) ID:", indiv, "]")

  per_y_ranints <- ranints[,indiv]
  y_pred_minus_ranints <- y_pred - per_y_ranints

  t(sapply(seq_len(nrow(y_pred_minus_ranints)),
           function(i) log(dnorm(y - y_pred_minus_ranints[i,]
                               , sd = sigma[i]))))
}

#' Get the fitted value from an ept model
#'
#' @param ept The ept model of interest
#' @return A n-vector of median fitted values
#' @export
fitted_ept <- function(ept) apply(extract(ept)$y_pred, 2, median)

#' Get the fitted values from a list of unpooled models
#'
#' @param np The data.frame of n no-pooling models
#' @return A n-vector of median fitted_values
#' @export
fitted_np <- function(np){
  lapply(np$fitted, function(m) fitted(m)) %>%
    unlist
}

#' Effect areas
#'
#' Create a posterior plot for the parameters with true effects
#' from results generated by `get_sglm_results` and `get_ept_results`
#'
#' @param results The result object of interest
#' @param effects The simulated effects
#' @param sds The node sds
#' @param tree The source tree
#' @param justLeaves Whether or not to just extract effects at leaves
#' @return a ggplot object
#' @export
effect_areas <- function(results, effects, sds, tree, justLeaves = FALSE){
  ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
  nodes <- tree$Get("name", filterFun = ff)
  
  scaled_effects <-
    effects[nodes] %>%
    `/`(sds[nodes]) %>%
    { data_frame(param = names(.), eff = .) }

  post <- results$b_post[,nodes] + results$fix[2]

  as <- mcmc_areas(post)

  as_data_u <-
    as$data %>%
    group_by(parameter) %>%
    slice(1) %>%
    ungroup %>%
    mutate(ymin = as.numeric(order(parameter, decreasing = TRUE))
         , ymax = ymin + .9
         , param = as.character(parameter)) %>%
    inner_join(scaled_effects, by = "param")

  as +
    geom_segment(aes(y = ymin, yend = ymax, x = eff, xend = eff)
               , col = "red", data = as_data_u) +
    scale_x_continuous()
}


#' Effect likelihoods
#'
#' Compute the likelihood of the simulated parameters given the models
#' 
#' @param results The result object of interest
#' @param effects The simulated effects
#' @param sds The node sds
#' @param tree The source tree
#' @param justLeaves just extract results for the leaves
#' @return A vector of pointwise likelihoods for each parameter
#' @export
pw_effect_loglik <- function(results, effects, sds, tree, justLeaves = FALSE){
  ff <- `if`(justLeaves, isLeaf, function(x) TRUE) 
  nodes <- tree$Get("name", filterFun = ff)
  effects <- effects[nodes]
  sds <- sds[nodes]
  
  scaled_effects <-
    effects %>%
    `/`(sds) %>%
    { data.frame(parameter = names(.), eff = .) }

  post <- results$b_post[,nodes] + results$fix[2]

  gauss_approx <-
    apply(post, 2, function(col) c(mean = mean(col), sd = sd(col)))

  imap_dbl(scaled_effects$eff
         , ~ log(.Machine$double.eps +
                 dnorm(.x, mean = gauss_approx[1, .y], sd = gauss_approx[2, .y])))
}

#' Compute the difference in beta between parent and child
#'
#' Add the posterior for the estimated effects back to a tree
#' then compute the difference between parent and child.
#'
#' @param results The result object of interest, must have b_post
#' and fix elements.
#' @param tree The hierarchy tree of interest
#' @return tree modified with extra node attributes b_post and b_diff
#' @export
compute_differences <-
  function(results, tree){
    ex_tree <- Clone(tree)

    ex_tree$Do(function(n) n$b_post <- results$b_post[,n$name])
    ex_tree$Do(function(n){
      if(is.null(n$parent)){
        n$b_diff <- n$b_post
      } else {
        n$b_diff <- n$b_post - n$parent$b_post
      }
    })
    
    ex_tree  
  }

