#' Modified Similarity Percentages (SIMPER) function
#'
#' Decomposing Bray-Curtis dissimilarity index to estimate the contribution of individual species within group
#'
#' @details \link{simper2} modifies from \link[vegan]{simper} function. The original \link[vegan]{simper} function only computes the discriminating species between two groups using Bray-Crutis dissimilarity.
#' The simper2 decomposes the within-group Bray-Crutis dissimilarity and computes the contribution of each species. One can also achieve the within-group SIMPER with the original \link[vegan]{simper} function by
#' copying the community data and creasting 2 dummy groups; however, the \link{simper2} should run faster. More details see \link[vegan]{simper}.
#' @param comm Community data matrix.
#' @param group Factor describing the group structure. Must have at least 2 levels.
#' @return A list of data frame with variables inculding:
#' \describe{
#' \item{overall}{The overall within-group dissimilarity}
#' \item{average}{Average contribution to overall dissimilarity}
#' \item{sd}{Standard deviation of contribution}
#' \item{ratio}{Average to sd ratio}
#' \item{avg}{Average abundances within group}
#' \item{contr}{\% contribution of speices}
#' \item{cusum}{Ordered cumulative contribution}
#' }
#' @export
simper2 <-
  function(comm, group = rep(1, nrow(comm))) {
    grp.no <- length(levels(factor(group)))
    if(grp.no > 1) { outlist <-list()
    p.no <- 0
    }
    for(i in levels(factor(group))){
      id <- which(group == i)
      comm1 <- comm[id, ]
      comp <- t(combn(1:nrow(comm1), 2))
      contr <- matrix(ncol = ncol(comm1), nrow = nrow(comp))
      for (j in 1:nrow(comp)) {
        md <- abs(comm1[comp[j, 1], ] - comm1[comp[j, 2], ])
        me <- comm1[comp[j, 1], ] + comm1[comp[j, 2], ]
        contr[j, ] <- md/sum(me)
      }
      average <- colMeans(contr)
      overall <- sum(average)
      sdi <- apply(contr, 2, sd)
      ratio <- average/sdi
      avg <- colMeans(comm1)
      ord <- order(average, decreasing = TRUE)
      cusum <- cumsum(average[ord]/overall)
      out <- cbind(average = average[ord], sd = sdi[ord],
                   avg = avg[ord], ratio = ratio[ord],
                   contr = average[ord]/overall, cusum = cusum)
      out <- list(overall = overall, summary = out)
      if(grp.no > 1) { p.no <- p.no+1
      outlist[[p.no]]<-out
      names(outlist)[[p.no]] <- i
      }
    }
    if(grp.no > 1) return(outlist) else return(out)
  }
