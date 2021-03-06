% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarization.R
\name{compute_differences}
\alias{compute_differences}
\title{Compute the difference in beta between parent and child}
\usage{
compute_differences(results, tree)
}
\arguments{
\item{results}{The result object of interest, must have b_post
and fix elements.}

\item{tree}{The hierarchy tree of interest}
}
\value{
tree modified with extra node attributes b_post and b_diff
}
\description{
Add the posterior for the estimated effects back to a tree
then compute the difference between parent and child.
}
