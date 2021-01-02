# use original fb networks

library(igraph)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(nett)

# compare different stats -------------------------------------------------
get_comm_profiles = function(A, 
                             get_labels = function(A, K) {
                               nett::spec_clust(A, K, tau = 0.25, niter = 20)
                               # nett::fastCPL(A, K, ilabels = spec_clust(A, K, tau=0.25))
                             },
                             Kmin = 1, Kmax = 25) {
  Ks <- Kmin:Kmax
  labels = sapply(Kmin:(Kmax+1), function(k) get_labels(A, k))
  
  tibble(
    "K" = Ks, 
    # "NAC" = sapply(Ks, function(k) nac_test(A, k, z=labels[ , k-Kmin+1], y = labels[, k-Kmin+2])$stat), 
    "SNAC" = sapply(Ks, function(k) snac_test(A, k, labels[ , k-Kmin+1])$stat), 
    # "AS (SBM)" = sapply(Ks, function(k) spec_test(A, k, labels[ , k-Kmin+1], DC = F)),
    # "AS" = sapply(Ks, function(k) spec_test(A, k, labels[ , k-Kmin+1], DC = T)), 
    # "BIC" = sapply(Ks, function(k) -eval_dcsbm_bic(A, z = labels[ , k-Kmin+1], k, poi = T))
  )
}

# Load FB data
if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))

# Compute and plot profiles -----------------------------------------------
library(parallel)
res_fname = file.path("FBhist","fb_profs_all.rds")
if (!file.exists(res_fname)) {
  prof_res = do.call(bind_rows, mclapply(1:100, function(net_id) {
    fname = file.path("FBhist", sprintf("fb_profs_%d.rds", net_id))
    if (file.exists(fname)) {
      temp = readRDS(fname)
    } else {
      A = Alist[[net_id]]    
      try({
        temp = get_comm_profiles(A) %>% mutate("net_id" = net_id)
        saveRDS(temp, fname)
      }, silent = T)
    }
    temp
  }, mc.cores = 20))
  saveRDS(prof_res, res_fname)
} else {
  prof_res = readRDS(res_fname)
}

# Generate DCSBM profles ---------------------------------------------------------
res_fname = file.path("FBhist","dcsbm_profs_all.rds")
Ktru = 3
Pi = rep(1,Ktru)/Ktru
if (!file.exists(res_fname)) {
  dt = system.time(
    dcsbm_prof_res <- do.call(bind_rows, mclapply(1:100, function(net_id) {
      fname = file.path("FBhist", sprintf("dcsbm_profs_%d.rds", net_id))
      if (file.exists(fname)) {
        temp = readRDS(fname)
      } else {
        n = nrow(Alist[[net_id]])
        ave_deg = sum(Alist[[net_id]])/n
        B <- pp_conn(n, 0.1, lambda = ave_deg, pri = Pi)$B
        theta <- EnvStats::rpareto(n, 3/4, 4)
        z <- sample(Ktru, n, replace=T, prob=Pi)
        A = nett::sample_dcsbm(z, B, theta = theta)
        try({
          temp = get_comm_profiles(A) %>% mutate("net_id" = net_id)
          saveRDS(temp, fname)
        }, silent = T)
      }
      temp
    }, mc.cores = 20))
  )["elapsed"]
  cat(sprintf("dt = %3.3f (s)", dt))
  saveRDS(dcsbm_prof_res, res_fname)
} else {
  dcsbm_prof_res = readRDS(res_fname)
}

both_prof_res = bind_rows( prof_res %>% add_column(data = "FB"), 
                           dcsbm_prof_res %>% add_column(data = "DCSBM") )

# The following is just to increase the vertical spacing between legend keys 
# See https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
GeomViolin$draw_key = function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

both_prof_res %>%  
  filter(K %in% c(1,2,3,4,5,10,20,15,25)) %>% 
  tidyr::pivot_wider(names_from = data, values_from = SNAC) %>% 
  ggplot(aes(x=factor(K))) + 
  geom_violin(aes(y=FB, color = "FB-100", fill="FB-100"), trim=T, alpha=0.5) +
  geom_violin(aes(y=DCSBM,  color = "DCSBM", fill="DCSBM"), trim=T, alpha=0.5) +
  # coord_cartesian(ylim=c(-5, 200)) + # do not use ylim(-5,200), causes data to be thrown out
  theme_minimal(base_size = 18) +
  xlab("Number of communities (K)") +
  ylab("Value of SNAC+ statistic") +
  scale_colour_manual(name="data", values=c("FB-100"="lightblue", DCSBM="red")) +
  scale_fill_manual(name="data", values=c("FB-100"="lightblue", DCSBM="#F8766D")) +
  theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.88, 0.9)
    # legend.text = ggplot2::element_text(size=25),
    # text = ggplot2::element_text(size=26)
  ) 

ggsave("FBhist/fb_hist.pdf",width=6, height=5)
ggsave("FBhist/fb_hist_zoom_in.pdf",width=6, height=5)




