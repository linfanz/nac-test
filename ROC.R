# ROC 
library(Matrix)
library(tibble)
library(dplyr)
library(ggplot2)
library(nett)

run_sims = T 
lvm_alt = T #whether to use DCLVM as the alternative model

nruns = 1000 # num of simulated adjacency matrices
DC <- T # whether the model has the degree-correction 
balanced <- F # whether community sizes are balanced
K <- 4 # number of communities
n <- 10000
oir <- 0.1 # out-in ratio
lambda <- 15 # average degrees
K_H0 <- K # number of communities in H0
K_H1 <- K+1 # number of communities in Ha

if (DC) {
  theta <- EnvStats::rpareto(n, 3/4, 4)    
} else {
  theta <- 1
}

if (balanced) {
  pri_func = function(K) rep(1, K)
  # pri0 = rep(1, K_H0)
  # pri1 = rep(1, K_H1)
} else {
  pri_func = function(K) 1:K
  # pri0 = 1:K_H0
  # pri1 = 1:K_H1
}


file_tag = gsub("\\.","p",sprintf("roc6_%dvs%d_n%d_lam%d_oir%2.2f_nruns%d_%s%s_%s", 
                                  K_H0, K_H1, n, lambda, oir, nruns, ifelse(DC, "dc", ""), ifelse(lvm_alt,"lvm","sbm"),  ifelse(balanced,"bal","unbal")))


apply_methods = function(A) {
  z0 = spec_clust(A, K_H0)
  z0p = spec_clust(A, K_H0+1)
  z1 = spec_clust(A, K_H1)
  tibble::tribble(
    ~method, ~tstat, ~twosided,
    "NAC+", nac_test(A, K_H0, z=z0, y=z0p)$stat, F,
    "SNAC+", snac_test(A, K_H0, z=z0)$stat, F,
    "AS", adj_spec_test(A, K_H0, z=z0), F,
    "LR", eval_dcsbm_loglr(A, cbind(z0, z1), poi=T), F
  )
}

gen_lvm = function(n, Ktru, lambda, theta) {
  d = Ktru
  labels = sample(Ktru, n, replace=T, prob= pri_func(Ktru))
  mu = diag(Ktru)
  z = 2*mu[labels, ] + 1*matrix(rnorm(n*d), n)
  sample_dclvm(z, lambda, theta)
}

gen_null_data = function() {
  sample_dcpp(n, lambda,  K_H0, oir = oir, theta, pri = pri_func(K_H0))$adj
  # gen_lvm(n, K_H0, lambda, theta)
}

gen_alt_data = function() {
  if (lvm_alt) {
    return( gen_lvm(n, K_H1, lambda, theta) )
  } else {
    return( sample_dcpp(n, lambda,  K_H1, oir = oir, theta, pri = pri_func(K_H1))$adj )
  }
}

if (run_sims) {
  roc_res = simulate_roc(apply_methods,
                         gen_null_data = gen_null_data,
                         gen_alt_data = gen_alt_data,
                         nruns = nruns, core_count = 24)
  printf('time = %3.3f', roc_res$elapsed_time)
  saveRDS(roc_res, paste0(file_tag, ".rds"))
} else {
  roc_res = readRDS(paste0(file_tag, ".rds"))
}

plot_roc <- function(roc_results, method_names=NULL) {
  if (!is.null(method_names)){
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method, levels = method_names))
  } else {
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method))
  }
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p = roc_results %>%
    ggplot2::ggplot(ggplot2::aes(x = FPR, y = TPR, color = method, linetype = method)) +
    ggplot2::scale_colour_manual(values=cbbPalette)+
    ggplot2::geom_line(size=2)  +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=18))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::geom_abline(intercept =0 , slope = 1, linetype="dashed") +
    ggplot2::scale_x_continuous(limits = c(0,1.01), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(0,1.01),expand = c(0,0)) +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.8, 0.2),
      legend.text = ggplot2::element_text(size=25),
      text = ggplot2::element_text(size=26)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 4, keyheight = 1)
    ) 
  p
}

plot_res <- roc_res$roc
plot_res$method <- factor(plot_res$method, levels = c("NAC+","SNAC+","AS", "LR"))
plot_roc(plot_res)
# ggsave(paste0(file_tag, ".pdf"), width = 8.22, height = 6.94)

