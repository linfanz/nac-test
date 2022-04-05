library(fastRG) # used to generate Poisson DCSBM
library(dplyr)
library(nett)
set.seed(123)

if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))


# modify snac_test to compute on large degree only   ----------------------
LargeDegIdx <- function(z, K, deg, perc=0.5) {
  return(
    as.vector(unlist(sapply(1:K, function(k) {
      index_k <- which(z == k)
      index_k[which(deg[index_k] > median(deg[index_k]))]
    })))
  )
}

snac_test_v2 <- function(A, K, z = NULL, ratio = 0.5, fromEachCommunity = TRUE, 
                      plus = TRUE, trueLabel = FALSE, useLargeDeg = FALSE,
                      cluster_fct = spec_clust, nrep = 1, ...) {
  if (is.null(z)) 
    z <- cluster_fct(A, K, ...)
  
  n <- length(z)
  stat <- c()
  deg <- rowSums(A)
  for (i in 1:nrep) {
    ## get index 
    if (fromEachCommunity) {
      if (useLargeDeg){
        index2 <- LargeDegIdx(z, K, deg)
      }else{
        index2 <- sampleEveryComm(z, K, ratio)
      }
    }
    else {
      if (useLargeDeg){
        index2 <- (1:n)[which(deg > median(deg))]
      }else{
        index2 <- sample(n, round(n * ratio))
      }
    }
    index1 <- (1:n)[-index2]
    
    ## get y label
    if (plus) {
      y1 <- cluster_fct(A[index1, index1], K + 1, ...)
    }
    else {
      if (trueLabel) {
        y1 <- z[index1]
      }else{
        y1 <- cluster_fct(A[index1, index1], K, ...)
      }
    }
    z2 <- z[index2]
    stat[i] <- nac_test(A[index2, index1], K, z = z2, y = y1)$stat
  }
  return(list(stat = stat, z = z))
}


# synthesize DCSBM  -------------------------------------------------------
Ktru <- 3
all_degs = unlist(lapply(Alist, rowSums))
# all_degs = all_degs[all_degs < 800]
N = length(all_degs)
# estimate pareto parameters 
location_h = min(all_degs)
shape_h = N/(sum(log(all_degs)) - N*log(location_h))

### kernel density estimation 
# bw <- density(all_degs)$bw

deg_match_res <- do.call(bind_rows, mclapply(1:100, function(net_id) {
  A_fb <- Alist[[net_id]]
  deg <- rowSums(A_fb)
  ave_deg <- mean(deg)
  n <- dim(A_fb)[1]
  
  ### a) use pareto distribution to generate theta
  # estimate pareto parameters for each individual network
  # location_h = min(deg)
  # shape_h = n/(sum(log(deg)) - n*log(location_h))
  ## or omit the above two lines and use all degrees to estimate parameter 
  # theta <- rpareto(n, location_h, shape_h)
  
  ### b) use kernel density estimation result
  # theta <- rnorm(n, mean = all_degs, sd = bw)
  # theta <- theta[theta > 0]
  # theta <- sample(theta, n, replace = TRUE)
  
  ### c) use FB degrees directly as theta
  theta = deg 
  
  B = pp_conn(n, oir = 0.1, lambda = ave_deg,
               pri = rep(1,Ktru), theta = theta)$B
  
  #Bernoulli generation
  z_ber <- sample(Ktru, n, replace = TRUE)
  A_ber <- nett::sample_dcsbm(z_ber, B, theta) 
  deg_ber <- rowSums(A_ber)
  # Poisson generation
  model <- fastRG::dcsbm(theta = theta, B = B, expected_degree = ave_deg)
  levels(model$z) <- 1:Ktru
  z_poi <- as.numeric(model$z)
  A_poi <- fastRG::sample_sparse(model)
  deg_poi <- rowSums(A_poi)
  
  # results 
  tibble::tribble(
    ~method, ~KS_pvalue, ~ave_deg_diff, ~SNAC_zh, ~SNAC_z, 
    "ber", ks.test(deg, deg_ber)$p.value, abs(mean(deg_ber) - ave_deg), 
    snac_test_v2(A_ber, Ktru, plus = F, useLargeDeg = T)$stat, 
    snac_test_v2(A_ber, Ktru, z = z_ber, plus = F, trueLabel = T, useLargeDeg = T)$stat,
    
    "poi", ks.test(deg, deg_poi)$p.value, abs(mean(deg_poi) - ave_deg),
    snac_test_v2(A_poi, Ktru, plus = F,  useLargeDeg = T)$stat, 
    snac_test_v2(A_poi, Ktru, z = z_poi, plus = F, trueLabel = T, useLargeDeg = T)$stat,
  )
}, mc.cores = 10))


# results -----------------------------------------------------------------

### synthesized DCSBM has degree distribution close to FB networks 
mean(deg_match_res$KS_pvalue > 0.01)

### empirical histograms
p1 = deg_match_res %>%
  filter(method == "ber")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC (true labels)")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p2 = deg_match_res %>%
  filter(method == "poi")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC (true labels)")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p3 = deg_match_res %>%
  filter(method == "ber")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC (est labels)")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p4 = deg_match_res %>%
  filter(method == "poi")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC (est labels)")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

ggpubr::ggarrange(p1, p2, p3, p4, 
                  labels = c("ber_z", "poi_z", "ber_zh", "poi_zh"),
                  ncol = 2, nrow = 2)

