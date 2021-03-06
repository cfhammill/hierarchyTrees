% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R
\name{h2}
\alias{h2}
\title{Fit a two-parent hierarchical model}
\usage{
h2(data, cores = 4)
}
\arguments{
\item{data}{A data object generated by \code{tree_to_edt_data}}

\item{cores}{number of cores to use, default 4}
}
\description{
Fits \code{scaled_vol ~ group + (group | p0) + (group | p1) + (group | p2) + (1 | ID)}
with \link{stan_lmer}
}
