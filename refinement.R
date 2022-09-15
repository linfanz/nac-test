# when y is refinement of z
nmi_super <- function(z,y){
  K = length(table(z))
  L = length(table(y))
  if (K >=  L){
    return ("Wrong!")
  }else{
    CM = t(label_vec2mat(z, K)) %*% label_vec2mat(y, L)
    normCM = CM/sum(CM)
    IDX = CM == 0
    jointEnt = -sum((normCM[!IDX]) * log(normCM[!IDX]))
    indpt = matrix(rowSums(normCM), ncol = 1) %*% matrix(colSums(normCM), nrow = 1)
    
    yCondz = normCM/matrix(rep(colSums(normCM), K), nrow =  K, byrow = T)
    zCondy = normCM/matrix(rep(rowSums(normCM), L), ncol = L)
    yzEnt = -sum(normCM[!IDX] * log(yCondz[!IDX]))
    zyEnt = -sum(normCM[!IDX] * log(zCondy[!IDX]))
    # jointEnt - yzEnt - zyEnt
    
    muInfo = sum(normCM[!IDX] * log(normCM[!IDX]/indpt[!IDX])) 
    
    # muInfo/(jointEnt - yzEnt)
    return(muInfo/(jointEnt - zyEnt))
  }
}

# muInfo/jointEnt
## test
# z = c(1,1,2,2,2,2)
# y = c(1,1,2,2,3,3)
# nmi_super(z,y)
# 
# n <- 10000 # number of nodes
# K <- 4 # number of communities
# beta <- 0.1 # out-in ratio
# lambda <- 30 # average degrees
# Pi <- rep(1,K)/K # distribution of cluster labels
# B <- pp_conn(n, beta, lambda, Pi)$B
# theta <- rpareto(n, 3/4, 4)
# 
# z <- sample(K, n, replace = T)
# A <- sample_dcsbm(z, B, theta = theta)
# zh <- spec_clust(A, K)
# nett::compute_mutual_info(z, zh)
# yh <- spec_clust(A, K+1)
# nmi_super(z, yh)
