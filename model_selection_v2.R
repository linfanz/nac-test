# model selection performance with sequential test
## add more cases 

library(nett)
library(EnvStats)
library(ggplot2)
library(Matrix)
library(dplyr)
library(parallel)
library(foreach)
library(randnet)

source('bootstrap debiasing.R')
 
# parameters set up -------------------------------------------------------
n <- 5000 # number of nodes
Ktru <- 4 # number of communities
Kmax <- 6 # maximum community number tested
Pi <- rep(1,Ktru)/Ktru # distribution of cluster labels
# Pi <- (1:Ktru)/sum(1:Ktru) # added unbalanced size v1
# Pi <- c(1,2,2,3)/8
# Pi <- c(1,1,2,3)/7
beta <- 0.3 # out-in-ratio for planted partition model
beta <- 0.2
lamax = 45
# average degree
lambdas <- round(10^seq(log10(5),log10(lamax), length.out = 10), 2)
nlam = length(lambdas)
# connectivity matrix 
# B1
Bs <- lapply(lambdas,  function(lambda) pp_conn(n, beta, lambda, Pi)$B )
# B2
# gamma <- 0.3
# B3,  different numbers on the diagonal
# Bs <- lapply(lambdas,  function(lambda) pp_conn(n, beta, lambda, Pi, d = c(1,2,3,1))$B)
# B4
# B4 <- matrix(c(0.9, 0.45, 0.5, 0.3,
#                  0.45, 0.6, 0.25, 0.2,
#                  0.5, 0.25, 0.1, 0.4,
#                0.3, 0.2, 0.4, 0.1), nrow = 4, byrow = T)
# Bs <- lapply(lambdas, function(lambda){
#   scale = get_dcsbm_exav_deg(n, Pi, B4, 1)
#   B4 * lambda/scale
# })
# Bs

# simulation --------------------------------------------------------------
bootrep <- 10 # bootstrap debaising repetition 
alpha <- 1e-6 # significance level
nitr <- 200 # number of simulated adjacency matrix
runs <- expand.grid(lam_idx = 1:nlam, itr = 1:nitr)

simulate_run <- function(j) {
  set.seed(j)
  idx <- runs[j,"lam_idx"]
  itr <- runs[j,"itr"]
  
  lambda <- lambdas[idx]
  
  # B1 or B3
  B <- Bs[[idx]]
  # B2
  # B = gen_rand_conn(n, Ktru, lambda, gamma)
  
  z <- sample(Ktru, n, replace=T, prob=Pi)
  theta <- rpareto(n, 3/4, 4)
  A <- nett::sample_dcsbm(z, B, theta = theta)
  
  # A <- nett::fast_sbm(z, B)
  # compute all the estimated labels first
  labels = sapply(1:(Kmax+1), function(k) nett::spec_clust(A, k))
  nac_results <- model_select(A, labels, alpha=alpha, bootrep=bootrep, Kmax=Kmax) 
  
  # all statistic
  res = tidyr::tribble(
    ~method, ~Kest,
    "SNAC+", nac_results$snac,
    "SNAC+ boot", nac_results$snac_boot,
    "NAC+ boot", nac_results$nac_boot,
    "AS boot", nac_results$as_boot,
    "BIC", nac_results$LRBIC,
    "BH",bethe_hessian_select(A, Kmax)$K, # BHMC
    "NCV", which.min(NCV.select(A,Kmax)$dc.l2), # network cv
    "ECV", which.min(ECV.block(A, Kmax)$dc.l2) # edge cv
  )
  res$lambda = lambda
  res$Ktru = Ktru
  res$iter = itr
  res$nmi = compute_mutual_info(z,  spec_clust(A, Ktru, tau=0.25)) 
  
  return(res)
}

CPU_CORES_TO_USE <- min(30, detectCores()-1)
printf('Using %d cores ...', CPU_CORES_TO_USE)

elapsed_time <- system.time( 
  result <- do.call(rbind, mclapply(1:nrow(runs), simulate_run, mc.cores = CPU_CORES_TO_USE))
)["elapsed"]
printf('Total time    = %3.1f (s)\nTime per iter = %3.1f (s)', elapsed_time, elapsed_time/nrow(runs))


# visualization -----------------------------------------------------------

nmi <- result %>% 
  group_by(lambda) %>% 
  summarize(nmi = mean(nmi))

acc_res <- result %>% 
  group_by(method,lambda) %>% 
  summarize(acc = mean(Kest==Ktru), acc_sd = sd(Kest==Ktru)/sqrt(nitr))

acc_res$method <- factor(acc_res$method,
                         levels = c("NAC+ boot","SNAC+ boot",  "AS boot", "SNAC+", "BIC", "BH", "ECV", "NCV"
                         ))
acc_res <- na.omit(acc_res)
acc_res$lambda <- as.numeric(acc_res$lambda)

## all methods
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
as_tibble(acc_res) %>% 
  filter(!grepl("th",method)) %>% 
  ggplot(aes(x=lambda, y=acc, col=method)) + 
  geom_line(aes(linetype=method), size=2) + 
  scale_linetype_manual(values=c(1,2,5,1,4,6,6,6))+
  scale_colour_manual(values=cbbPalette)+
  theme_bw()+
  geom_point(aes(shape=method), size = 5)+
  scale_shape_manual(
    values = c(15, 17, 23, 13, 19, 25, 18, 16)
  ) +
  # geom_errorbar(aes(ymin=acc-acc_sd, ymax=acc+acc_sd), width=.2)+
  ylab("Expected Accuracy") + xlab("Average Degree")+
  theme(legend.position="none", text = element_text(size=26)) # use this to remove legend
  theme(legend.background=element_blank(),
        legend.title=element_blank(),
        legend.position = c(0.83, 0.5),
        legend.text = element_text(size=25),
        text = element_text(size=26),
        legend.key.width = unit(4, "line")) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch"))


# save results
base_fname = gsub("\\.","p", sprintf("seq_BerDCSBM_Ktru%d_Kmax%d_n%d_beta%2.2f_m%d_B1_unbalanced_v1", Ktru, Kmax, n, beta, nitr))
ggsave(paste0(base_fname,".pdf"), width=9.92, height=8.15)
# saveRDS(result, paste0("result_", base_fname,".rds"))
  
# base_fname = gsub("\\.","p", sprintf("seq_BerDCSBM_Ktru%d_Kmax%d_n%d_beta%2.2f_m%d_B1", Ktru, Kmax, n, beta, nitr))
# ggsave(paste0(base_fname,".pdf"), width=9.92, height=8.15)
# saveRDS(result, paste0("result_", base_fname,".rds"))

seq_BerDCSBM_Ktru4_Kmax6_n5000_beta0p30_m200_B1