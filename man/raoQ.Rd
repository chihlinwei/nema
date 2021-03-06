% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raoQ.R
\name{raoQ}
\alias{raoQ}
\title{Rao's quadratic entropy Q for pooled assemblages of N sites}
\source{
\url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}
}
\usage{
raoQ(data, dij, datatype)
}
\arguments{
\item{data}{a list with N sites; each element of list is species abundances or species-by-sampling-unit incidence matrix.}

\item{dij}{a matrix of species-pairwise distances.}

\item{datatype}{data type of input data: individual-based abundance data (datatype = "abundance") or species-by-sampling-unit incidence matrix (datatype = "incidence_raw").}
}
\value{
a numeric number of Q for pooled assemblages.
}
\description{
Calculating Rao's quadratic entropy Q for pooled assemblages of N sites based on Chao et al. (2018)
}
\details{
raoQ is a function of obtaining Rao's quadratic entropy Q. The original R code is ported from \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}.
The function is implemented here for convenience. The details should refer to Chao et al. (2018).
}
\references{
Chao et al. (2018) An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures.
Under revision, Ecological Monographs.
}
