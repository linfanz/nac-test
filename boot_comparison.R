# compare bootstrap methods 

# model selection performance with sequential test
library(nett)
library(EnvStats)
library(ggplot2)
library(Matrix)
library(dplyr)
library(parallel)
library(foreach)
library(randnet)

source('bootstrap_all_methods.R')

# parameters set up -------------------------------------------------------
n <- 5000 # number of nodes
Ktru <- 4 # number of communities
Kmax <- 6 # maximum community number tested
Pi <- rep(1,Ktru)/Ktru # distribution of cluster labels
#Pi <- c(1,1,2,3)/7
beta <- 0.2 # out-in-ratio for planted partition model 
lamax = 25 
# average degree
lambdas <- round(10^seq(log10(3),log10(lamax), length.out = 10), 2)
nlam = length(lambdas)
# connectivity matrix 
# B1
Bs <- lapply(lambdas,  function(lambda) pp_conn(n, beta, lambda, Pi)$B )
# B2
# gamma <- 0.3
# B3,  different numbers on the diagonal
# Bs <- lapply(lambdas,  function(lambda) pp_conn(n, beta, lambda, Pi, d = c(1,2,3,1))$B)

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
  
  theta <- rpareto(n, 3/4, 4)
  z <- sample(Ktru, n, replace=T, prob=Pi)
  
  A <- nett::sample_dcsbm(z, B, theta = theta)
  # compute all the estimated labels first
  labels = sapply(1:(Kmax+1), function(k) nett::spec_clust(A, k))
  
  nac_results <- nac_model_select(A, labels, alpha=alpha, bootrep=bootrep, Kmax=Kmax) # NAC methods
  
  # all statistic
  res = tidyr::tribble(
    ~method, ~Kest,
    "SNAC+", nac_results$snac,
    "SNAC+ boot de", nac_results$snac_boot_de,
    "FNAC+ boot de", nac_results$nac_boot_de,
    "SNAC+ boot qu", nac_results$snac_boot_qu,
    "FNAC+ boot qu", nac_results$nac_boot_qu,
    "SNAC+ dc boot de", nac_results$snac_dc_boot_de,
    "FNAC+ dc boot de", nac_results$nac_dc_boot_de,
    "SNAC+ dc boot qu", nac_results$snac_dc_boot_qu,
    "FNAC+ dc boot qu", nac_results$nac_dc_boot_qu
  )
  res$lambda = lambda
  res$Ktru = Ktru
  res$iter = itr
  res$nmi = compute_mutual_info(z,  labels[ , Ktru]) 
  
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
                         levels = c(    "SNAC+", 
                                        "SNAC+ boot de", 
                                        "FNAC+ boot de", 
                                        "SNAC+ boot qu",
                                        "FNAC+ boot qu", 
                                        "SNAC+ dc boot de", 
                                        "FNAC+ dc boot de", 
                                        "SNAC+ dc boot qu", 
                                        "FNAC+ dc boot qu"
                         ))
acc_res$lambda <- as.numeric(acc_res$lambda)

## all methods
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red")

# NAC version
as_tibble(acc_res) %>% 
  filter(!grepl("th",method) & method %in% c("SNAC+",
                                              "FNAC+ boot de", 
                                             "FNAC+ boot qu", 
                                             "FNAC+ dc boot de", 
                                             "FNAC+ dc boot qu")) %>% 
  ggplot(aes(x=lambda, y=acc, col=method)) + 
  geom_line(aes(linetype=method), size=2) + 
  scale_linetype_manual(values=c(1,2,5,4,6))+
  scale_colour_manual(values=cbbPalette[1:5])+
  theme_bw()+
  geom_point(aes(shape=method), size = 5)+
  scale_shape_manual(
    values = c(15, 17, 23, 13, 19)
  )+
  theme(text = element_text(size=20))+
  # geom_errorbar(aes(ymin=acc-acc_sd, ymax=acc+acc_sd), width=.2)+
  ylab("Expected Accuracy") + xlab("Average Degree")+
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.75, 0.3),
        legend.text = element_text(size=25),
        text = element_text(size=25),
        legend.key.width = unit(4, "line")) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch"))

# save results
# ggsave("nac_boot_comparison_1_z.pdf")
ggsave(paste0("fnac_boot_comparison_1_zh_(1).pdf"), width=9.92, height=8.15)

# snac version
as_tibble(acc_res) %>% 
  filter(!grepl("th",method) & method %in% c("SNAC+", "SNAC+ boot de", "SNAC+ boot qu",
                                             "SNAC+ dc boot de",
                                             "SNAC+ dc boot qu")) %>%
  ggplot(aes(x=lambda, y=acc, col=method)) + 
  geom_line(aes(linetype=method), size=2) + 
  scale_linetype_manual(values=c(1,2,5,1,4,6,6,6,6))+
  scale_colour_manual(values=cbbPalette[c(1, 6:9)])+
  theme_bw()+
  geom_point(aes(shape=method), size = 5)+
  scale_shape_manual(
    values = c(15, 25, 18, 16, 12)
  )+
  theme(text = element_text(size=20))+
  # geom_errorbar(aes(ymin=acc-acc_sd, ymax=acc+acc_sd), width=.2)+
  ylab("Expected Accuracy") + xlab("Average Degree")+
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.75, 0.3),
        legend.text = element_text(size=25),
        text = element_text(size=25),
        legend.key.width = unit(4, "line")) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch"))

# save results
# ggsave("snac_boot_comparison_1_z.pdf")
ggsave(paste0("snac_boot_comparison_1_zh_(1).pdf"), width=9.92, height=8.15)

# base_fname = gsub("\\.","p", sprintf("seq_BerDCSBM_Ktru%d_Kmax%d_n%d_beta%2.2f_m%d", Ktru, Kmax, n, beta, nitr))
# ggsave(paste0(base_fname,".pdf"))
# saveRDS(result, paste0("result_", base_fname,".rds"))