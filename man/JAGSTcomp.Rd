% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jagst.R
\name{JAGSTcomp}
\alias{JAGSTcomp}
\title{JAGST competitive test}
\usage{
JAGSTcomp(array, design, inds, its = 200, num_nest = 2000)
}
\arguments{
\item{array}{expression array. Genes along rows.}

\item{design}{design matrix (intercept included by default)}

\item{inds}{indices of transcripts to be tested}

\item{its}{number of iterations for null generation.  This number times num_nest is the total number of samples for the null.}

\item{num_nest}{number of samples within each null iteration.  Cheap computationally, but decreasing returns efficiency gain.}
}
\description{
JAGST competitive test
}
\examples{
arr <- matrix(rnorm(100*100),nrow=100,ncol=100)
des <- rep(0:1,each=50)
ind <- sample(100,10)
JAGSTcomp(arr,des,ind)
}
