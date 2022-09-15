# empirical distribution of SNAC and SNAC+ 

library(ggplot2)
library(EnvStats)
library(nett)
library(dplyr)
library(fastRG) # used to simulate Poisson DCSBM

n <- 10000 # number of nodes
K <- 4 # number of communities
beta <- 0.1 # out-in ratio
lambda <- 30 # average degrees
Pi <- rep(1,K)/K # distribution of cluster labels
B <- pp_conn(n, beta, lambda, Pi)$B
m <- 1000 # number of simulated adjecency matrix
z <- sample(K, n, replace=T, prob=Pi)
theta <- rpareto(n, 3/4, 4)
A <- sample_dcsbm(z, B, theta = theta)
zh <- spec_clust(A, K)
compute_mutual_info(z, zh)

boot_nac <- function(A, z, K, nrep = 500, dcsbm_boot = T) {
  # printf('nboot = %d ', nrep)
  # bootstrap with SBM
  out <- nett::estim_dcsbm(A, z)
  B <- out$B
  theta <- out$theta
  n <- length(theta)
  ave_deg <- sum(A)/n
  # bootstrap with DCSBM
  # theta <- out$theta
  boot_samp <- sapply(1:nrep, function(i) {
    # Aboot <- nett::fast_sbm(z, B)
    model_sbm <- sbm(n, K, B, pi = table(z)/n, expected_degree = ave_deg)
    # zboot <- spec_clust(Aboot, K)
    Aboot <- sample_sparse(model_sbm)
    
    model_dcsbm <- dcsbm(theta = theta, B = B, pi = table(z)/n, expected_degree = ave_deg)
    Aboot_DC <- sample_sparse(model_dcsbm)
    # Aboot_DC <- nett::sample_dcsbm(z, B, theta = theta)
    # zboot_DC <- spec_clust(Aboot_DC, K)
    c(
      snac_test(Aboot, K)$stat,
      snac_test(Aboot_DC, K)$stat
      # adj_spec_test(Aboot, K, zboot)
    )
  })
  return(boot_samp)
}

res <- boot_nac(A, zh, K)
res <- data.frame(SBM = res[1, ],  DCSBM = res[2, ])

ggplot(res, aes(x = SBM)) + 
  theme_bw() +  xlab("SNAC+")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white",
                 bins = 20)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )
# ggsave(paste0("SBM_boot_hist.pdf"), width = 8.22, height = 6.94)

ggplot(res, aes(x = DCSBM)) + 
  theme_bw() +  xlab("SNAC+")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white",
                 bins = 20)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )
# ggsave(paste0("DCSBM_boot_hist.pdf"), width = 8.22, height = 6.94)
