# empirical distribution of SNAC and SNAC+ 

library(ggplot2)
library(EnvStats)
library(nett)
library(dplyr)
library(fastRG) # used to simulate Poisson DCSBM

set.seed(123)
n <- 10000 # number of nodes
K <- 3 # number of communities
beta <- 0.05 # out-in ratio
lambda <- 40 # average degrees
Pi <- rep(1,K)/K # distribution of cluster labels
B <- pp_conn(n, beta, lambda, Pi)$B
z <- sample(K, n, replace=T, prob=Pi)
theta <- rpareto(n, 3/4, 4)

n <- length(theta)
ave_deg <- lambda 

boot_nac <- function(z, K, nrep = 2000, dcsbm_boot = T) {
  boot_samp <- sapply(1:nrep, function(i) {
    Aboot_ber <- nett::fast_sbm(z, B)
    model_sbm <- sbm(n, K, B, expected_degree = ave_deg)
    # # zboot <- spec_clust(Aboot, K)
    Aboot <- sample_sparse(model_sbm)
    
    model_dcsbm <- dcsbm(theta = theta, B = B, expected_degree = ave_deg)
    Aboot_DC <- sample_sparse(model_dcsbm)
    Aboot_DC_ber <- nett::sample_dcsbm(z, B, theta = theta)
    # zboot_DC <- spec_clust(Aboot_DC, K)
    c(
      snac_test(Aboot, K)$stat,
      snac_test(Aboot_ber, K)$stat,
      snac_test(Aboot_DC, K)$stat,
      snac_test(Aboot_DC_ber, K)$stat
      # adj_spec_test(Aboot, K, zboot)
    )
  })
  return(boot_samp)
}

res <- boot_nac(z, K)
res <- data.frame(SBM = res[1, ], berSBM = res[2, ],  DCSBM = res[3, ], berDCSBM = res[4, ] )

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
# ggsave(paste0("SBM_null_hist.pdf"), width = 8.22, height = 6.94)

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
# ggsave(paste0("DCSBM_null_hist.pdf"), width = 8.22, height = 6.94)

ggplot(res, aes(x = berSBM)) + 
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
# ggsave(paste0("berSBM_null_hist.pdf"), width = 8.22, height = 6.94)

ggplot(res, aes(x = berDCSBM)) + 
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
# ggsave(paste0("berDCSBM_null_hist.pdf"), width = 8.22, height = 6.94)


mean(res$DCSBM)
#-0.0477519
sd(res$DCSBM)
# 1.043478

mean(res$SBM)
# [1] -0.08605762
sd(res$SBM)
# [1] 0.9863768

shapiro.test(res$DCSBM)
shapiro.test(res$SBM)

ks.test(res$DCSBM, "pnorm")
#p-value = 0.01563
ks.test(res$SBM, "pnorm")
# p-value = 0.009916
## strangely, the KS test has large pvalue when there are 
## less data points  


# Bernoulli version
ks.test(res$berDCSBM, "pnorm")
ks.test(res$berSBM, "pnorm")
shapiro.test(res$berDCSBM)
shapiro.test(res$berSBM)

mean(res$berDCSBM)
# [1] -0.2915549
sd(res$berDCSBM)
# [1] 1.022245
mean(res$berSBM)
# [1] -0.2571552
sd(res$berSBM)
# [1] 0.9948815
