#' Rao's quadratic entropy Q for pooled assemblages of N sites
#'
#' Calculating Rao's quadratic entropy Q for pooled assemblages of N sites based on Chao et al. (2018)
#'
#' @details raoQ is a function of obtaining Rao's quadratic entropy Q. The original R code is ported from \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}.
#' The function is implemented here for convenience. The details should refer to Chao et al. (2018).
#'
#' @param data a list with N sites; each element of list is species abundances or species-by-sampling-unit incidence matrix.
#' @param dij a matrix of species-pairwise distances.
#' @param datatype data type of input data: individual-based abundance data (datatype = "abundance") or species-by-sampling-unit incidence matrix (datatype = "incidence_raw").
#' @return a numeric number of Q for pooled assemblages.
#' @references Chao et al. (2018) An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures.
#' Under revision, Ecological Monographs.
#' @source \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}
#' @export
raoQ <- function(data, dij, datatype){
  if (datatype=="abundance") {
    if(length(data)!=1){
      tmp <- apply(do.call(cbind,lapply(data, FUN = function(x) x/sum(x))), 1, mean)
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)
    }else{
      tmp <- data[[1]]/sum(data[[1]])
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)
    }
  } else {
    if(length(data)!=1){
      tmp <- apply(do.call(cbind,lapply(data, FUN = function(x) x[-1]/sum(x[-1]))), 1, mean)
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)
    }else{
      tmp <- data[[1]][-1]/sum(data[[1]][-1])
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)
    }
  }
  return(dmean)
}
