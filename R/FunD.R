#' Functional Diversity of N sites for various values of tau and q
#'
#' Calculating functional diversity of N sites for various values of tau and q based on Chao et al. (2018)
#'
#' @details FunD is a function of obtaining FD tau-table or q-table. The original R code is ported from \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}.
#' The function is implemented here for convenience. The details should refer to Chao et al. (2018).
#'
#' @param data a list with N sites; each element of list is species abundances or species-by-sampling-unit incidence matrix.
#' @param dij a matrix of species-pairwise distances.
#' @param tau a numeric or a vector of levels of threshold distinctiveness.
#' @param q a numeric or a vector of diversity orders; the suggested range for q is [0, 3].
#' @param boot a numeric of number of bootstrap replications.
#' @param datatype data type of input data: individual-based abundance data (datatype = "abundance") or species-by-sampling-unit incidence matrix (datatype = "incidence_raw").
#' @return two matrices of FD; the first matrix is the q-profile information, the second matrix is the tau profile information.
#' @references Chao et al. (2018) An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures.
#' Under revision, Ecological Monographs.
#' @source \url{https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt}
#' @export
FunD <- function(data, dij, tau, q, boot, datatype){
  EstiBootComm.Func_abun = function(data, distance, datatype){
    if (datatype=="abundance") {
      n = sum(data)
    } else if (datatype=="incidence_raw") {
      n <- data[1]
      data <- data[-1]
      u=sum(data)
    }
    distance = as.matrix(distance)
    dij = distance[data!=0, data!=0]

    X = data[data>0]
    f1 <- sum(X == 1) ; f2 <- sum(X == 2)
    f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
    if (datatype=="abundance") {
      C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
      W <- (1 - C1)/sum(X/n*(1-X/n)^n)
      Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)
    } else if (datatype=="incidence_raw") {
      C1 = ifelse(f2>0, 1-f1/u*(n-1)*f1/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/u/((n-1)*(f1-1)+2))
      W <- (1 - C1)/sum(X/u*(1-X/n)^n)
      Prob.hat.Unse <- rep(u/n*(1-C1)/f0.hat, f0.hat)
    }

    Prob.hat <- X/n*(1-W*(1-X/n)^n)
    Prob <- c(Prob.hat, Prob.hat.Unse)

    F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
    F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
    #
    if (datatype=="abundance") {
      F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
    } else if (datatype=="incidence_raw") {
      F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
    }
    if (f0.hat==0) {
      d=dij
    } else if (f0.hat==1) {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
      d00 = matrix(0, f0.hat, f0.hat)
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    } else {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)

      fo.num = (f0.hat * (f0.hat-1) )/2
      random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
      d00 = matrix(0, f0.hat, f0.hat)
      d00[upper.tri(d00)] = (F00hat/2)*random_d00
      d00 <- pmax(d00, t(d00))###signmatrix
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    }
    return(list("pi" = Prob,"dij" = d))
  }
  dij <-  as.matrix(dij)
  out <- as.vector(dij)
  out <- out[out!=0]
  dmin <- min(out)
  dmax <- max(out)
  if (datatype=="incidence_raw") {
    data <- lapply(data, function(i) {
      c(ncol(i), rowSums(i))
    })
  }
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
  FD.CI = function(data, dij, tau, q, datatype){
    if (datatype == "abundance") {
      qFun = FD_MLE(data, dij, tau, q)
    } else {
      qFun = FD_MLE(data[-1], dij, tau, q)
    }
    if(boot!=0){
      BT = EstiBootComm.Func_abun(data, dij, datatype)
      p_hat = BT[[1]]
      dij_boot = BT[[2]]
      dij_boot <-  as.matrix(dij_boot)
      dij_boot <- replace(dij_boot, dij_boot==0, 10^(-10))
      for (i in seq_len(nrow(dij_boot))) {
        dij_boot[i, i] <- 0
      }
      if (datatype=="abundance") {
        n=sum(data)
        Boot.X = rmultinom(boot, n, p_hat)
      } else {
        n=data[1]
        Boot.X = t(sapply(p_hat,function(i) rbinom(boot, n, i)))
      }
      qFun_sd = sd(sapply(seq_len(ncol(Boot.X)), function(i) {
        FD_MLE(Boot.X[, i], dij_boot, tau, q)
      }))
    }else{
      qFun_sd = 0
    }
    LCL = max(0, qFun - qnorm(0.975) * qFun_sd)
    UCL = qFun + qnorm(0.975) * qFun_sd
    a = round(c(qFun, qFun_sd, LCL, UCL), 4)
    a
  }

  Funq <- function(data, datatype){
    dminFDforq <- t(sapply(q, FUN = function(q) FD.CI(data, dij, dmin, q, datatype) ))
    dmaxFDforq <- t(sapply(q, FUN = function(q) FD.CI(data, dij, dmax, q, datatype) ))
    dmeanFDforq <-t(sapply(q, FUN = function(q) FD.CI(data, dij, dmean, q, datatype) ))
    out <- data.frame(rep(q,3), rbind(dminFDforq,dmaxFDforq,dmeanFDforq),rep(c("dmin","dmax","dmean"),each=length(q)))
  }
  Funtau <- function(data, datatype){
    q0FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 0, datatype) ))
    q1FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 1, datatype) ))
    q2FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 2, datatype) ))
    out <- data.frame(rep(tau,3), rbind(q0FDfortau,q1FDfortau,q2FDfortau),rep(c("0","1","2"),each=length(tau)))
  }

  if(length(data)!=1){
    name = names(data)
    Outputforq <- data.frame(do.call(rbind,lapply(data, Funq, datatype=datatype)), rep(name, each=3*length(q)), row.names = NULL)
    Outputfortau <- data.frame(do.call(rbind,lapply(data, Funtau, datatype=datatype)), rep(name, each=3*length(tau)), row.names = NULL)
  }else{
    name = names(data)
    Outputforq <- data.frame(Funq(data[[1]], datatype), name, row.names = NULL)
    Outputfortau <- data.frame(Funtau(data[[1]], datatype), name, row.names = NULL)
  }
  colnames(Outputforq) <- c("q","estimate", "s.e.", "LCL", "UCL", "tau","site")
  colnames(Outputfortau) <- c("tau","estimate", "s.e.", "LCL", "UCL", "q","site")

  Output <- list(forq = Outputforq, fortau = Outputfortau)
  return(Output)
}
