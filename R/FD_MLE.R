#' Functional Diversity of a single site for specified values of tau and q
#'
#' Calculating functional diversity of a single site for specified values of tau and q based on Chao et al. (2018)
#'
#' @details FD_MLE is a function of obtaining FD index of order q. The original R code is ported from \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}.
#' The function is implemented here for convenience. The details should refer to Chao et al. (2018).
#'
#' @param data a vector of species sample frequencies.
#' @param dij a matrix of species-pairwise distances.
#' @param tau a numeric for a specified level of threshold distinctiveness.
#' @param q a numeric for a specified diversity order q.
#' @return a numeric value of FD.
#' @references Chao et al. (2018) An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures.
#' Under revision, Ecological Monographs.
#' @source \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}
#' @export
FD_MLE <- function(data, dij, tau, q){
  dij <- as.matrix(dij)
  dij[which(dij>tau,arr.ind = T)] <- tau
  a <- as.vector((1 - dij/tau) %*% data )
  data <- data[a!=0]
  a <- a[a!=0]
  v <- data/a
  nplus <- sum(data)
  if(q==1){
    exp(sum(-v*a/nplus*log(a/nplus)))
  }else{
    (sum(v*(a/nplus)^q))^(1 / (1-q))
  }
}
