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
  function(model = c("effect_diffusion_tree", "simple_hierarchical")){
    model <- match.arg(model)
    stan_model(
      system.file(
        file.path("models",
                  paste0(model, ".stan"))
      , package = "hierarchyTrees")
    )
  }


#' Check the data requirements for a stan model
#'
#' Parse out the data block from a stan model and
#' report the required data
#'
#' @param model The stan model to inspect
#' @return a character vector of required data
#' @export
get_data_requirements <-
  function(model){
    split1 <- function(...) unlist(strsplit(...))

    ## Extract code Remove Comments
    code <-
      model@model_code %>%
      split1("\n") %>%
      sub("/.*", "", .) %>%
      paste0(collapse = "\n") 

    ## Find the data block
    data_block <-
      code %>%
      sub(".*?data\\{", "", ., perl = TRUE) %>%
      split1("") %>%
      Reduce(function(acc, ch){
        if(acc$done)    # If the data block is done return
          return(acc)   # sadly executes for each remaining char
        
        if(ch == "{")   # Notice if an internal block opens
          acc$open <- acc$open + 1

        if(ch == "}"){
          if(acc$open > 0){ # Notice if an internal block closes
            acc$open <- acc$open - 1
          } else {          # Notice the end of the data block
            acc$done <- TRUE
          }
        } else {            # Accumulate a data block character
          acc$stack <- c(acc$stack, ch)
        }

                           # return the accumulator
        return(acc)       
      }
    , .                         
    , list(stack = character(0) # empty char stack
         , open = 0             # 0 internal blocks open
         , done = FALSE))       # not done yet


    ## Reassemble the data block
    data_block <-
      data_block$stack %>%
      paste0(collapse = "") %>%
      split1("\n") %>%
      sub("^ +", "", .) %>%
      sub(" *; *", "", .) %>%
      Filter(function(st) st != "", .) 

    sub("([^ ]*) ([^[]*)(.*)", "\\1|\\2|\\3", data_block) %>%
      strsplit("\\|") %>%
      lapply(function(strs)
        data.frame(type = strs[1], name = strs[2], arr_size = strs[3]
                   , stringsAsFactors = FALSE)
        ) %>%
      bind_rows
  }

#' Create a data skeleton for a model
#'
#' Take the results of [get_data_requirements]
#' and produce a skeleton list object to be passed to
#' stan. If the requirements table has a value column
#' it will be include as the contents of the list. Otherwise
#' list elements are initialized to NA.
#'
#' @param req The data requirements table
#' @return A list with length equal to the number of rows
#' in `reqs` with names equal to the `name` column of
#' `reqs` and values initialized to either the contents
#' of the `value` column if it exists or NAs.
#' @md
#' @export
requirements_to_skeleton <-
  function(reqs){
    if(is.null(reqs$value))
      reqs$value <- NA

    reqs$value %>%
      as.list %>%
      setNames(reqs$name)
  }


#' Fix node names and plot a tree
#'
#' Fix node names then plot the tree
#' 
#' @param tree A tree to plot
#' @return NULL invisibly
#' @export
fix_names_and_plot <- function(tree){
  ex_tree <- Clone(tree)
  nodes <- ex_tree$Get("name")
  new_names <- setNames(make.names(nodes, unique = TRUE), nodes)

  fix_node_names <- function(tree, new_names){
    tree$name <- new_names[tree$name]
    if(is.null(tree$children) || length(tree$children) == 0)
      return(invisible(NULL))

    walk(tree$children, fix_node_names, new_names)
  }

  fix_node_names(ex_tree, new_names)

  plot(ex_tree, output = "visNetwork")
}
