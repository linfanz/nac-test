# quantile filtering version 
quantile_filter <- function(z, idx, deg, K, sigma) {
  return(
    as.vector(unlist(sapply(1:K, function(k) {
      pos_k <- which(z == k)
      th_k <- quantile(deg[pos_k], sigma)
      index_k <- idx[pos_k]
      index_k[which(deg[pos_k] > th_k)]
    })))
  )
}

sampleEveryComm <- function(z, K, ratio=0.5) {
  return(
    as.vector( unlist(sapply(1:K, function(k) {
      index_k <- which(z == k)
      sample(index_k, size = round(length(index_k)*ratio))
    })))
  )
}

snac_test_v2 <- function (A, K, z = NULL, ratio = 0.5, 
                          sigma = 0, 
                          fromEachCommunity = TRUE, plus = TRUE, 
                          cluster_fct = spec_clust, nrep = 1, ...) 
{
  if (is.null(z)) 
    z <- cluster_fct(A, K, ...)
  n <- length(z)
  stat <- c()
  for (i in 1:nrep) {
    if (fromEachCommunity) {
      index1 <- sampleEveryComm(z, K, ratio)
    }else{
      index1 <- sample(n, round(n * ratio))
    }
    index2 <- (1:n)[-index1]
    if (plus) {
      y1 <- cluster_fct(A[index1, index1], K + 1, ...)
    }else{
      y1 <- cluster_fct(A[index1, index1], K, ...)
    }
    z2 <- z[index2]
    if (sigma){
      degs2 <- rowSums(A[index2, index1])
      index2 <- quantile_filter(z2, index2, degs2, K, sigma)
      z2 <- z[index2]
    }
    
    stat[i] <- nac_test(A[index2, index1], K, z = z2, y = y1)$stat
  }
  return(list(stat = stat, z = z))
}
