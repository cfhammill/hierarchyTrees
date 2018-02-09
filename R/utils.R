#' Compile a tree model
#'
#' Compile one of the models available in this
#' package.
#' @param model The model to compile, see the
#' default args for a list of available models.
#' @return A [stan_model] object
#' @md
#' @export
compile_models <-
  function(model = c("effect_propagation_tree")){
    model <- match.arg(model)
    stan_model(
      system.file(
        file.path("models",
                  paste(model, ".stan"))
      , package = "hierarchyTreeTests")
    )
  }
