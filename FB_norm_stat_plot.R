library(igraph)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(nett)
library(tidyverse)
library(foreach)
library(doRNG)

compute_all_stats = function(A, 
                             get_labels = function(A, K) {
                               spec_clust(A, K, tau=0.25, niter = 20)
                               #nett::fastCPL(A, K, ilabels = spec_clust(A, K, tau=0.25))
                             },
                             Kmin = 1, Kmax = 13) {
  Ks <- Kmin:Kmax
  labels = sapply(Kmin:(Kmax+1), function(k) get_labels(A, k))
  list(
    labels = labels, 
    profiles = tibble(
      "K" = Ks, 
      "NAC" = sapply(Ks, function(k) nac_test(A, k, z=labels[ , k-Kmin+1], y = labels[, k-Kmin+2])$stat), 
      "SNAC" = sapply(Ks, function(k) snac_test(A, k, z = labels[ , k-Kmin+1])$stat), 
      "AS (SBM)" = sapply(Ks, function(k) adj_spec_test(A, k, labels[ , k-Kmin+1], DC = F)),
      "AS" = sapply(Ks, function(k) adj_spec_test(A, k, labels[ , k-Kmin+1], DC = T)), 
      "- BIC" = sapply(Ks, function(k) -eval_dcsbm_bic(A, z = labels[ , k-Kmin+1], k, poi = T))
    )
  )
}

# explore network -------------------------------------------------

if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))
college_names <- names(Alist)
# net_id = which(college_names == "Wellesley22")
net_id = which(college_names == "Harvard1")
# net_id = which(college_names == "Bucknell39")
g = graph_from_adjacency_matrix(Alist[[net_id]], mode = "undirected")
g = extract_largest_cc(g)
net_name = college_names[net_id]

degs = degree(g)
deg_prec = 0.75
g2 = g %>% 
  induced_subgraph(degs <= quantile(degs, deg_prec))%>%
  extract_largest_cc()
degs2 = degree(g2)
summary(degs)
summary(degs2)
plot_deg_dist(g2) # plot degree distribution
printf('# of nodes = %d
reduction in # of nodes = %2.1fx
reduction in max deg. = %dx
reduction in mean deg. = %dx', 
       vcount(g2), 
       vcount(g)/vcount(g2), 
       round(max(degs) / max(degs2)),
       round(mean(degs) / mean(degs2)))

A = as_adj(g2)

# Compute and plot normalized stats -----------------------------------------------
out2 = compute_all_stats(A)
profs = out2$profiles
# scale statistics to (-1,1) scale
profs_n = profs %>% 
  mutate_at(vars(-K), ~( . / max(abs(.))) ) 
profs_nl = profs_n %>% 
  pivot_longer(-K, names_to = "method", values_to = "Normalized Stat.")  %>% 
  mutate(method = replace(method, method == "SNAC", "SNAC+")) %>% 
  mutate(method = replace(method, method == "NAC", "NAC+")) %>%
  mutate(method = factor(method, c("SNAC+","NAC+","- BIC","AS", "AS (SBM)")))

profs_nl %>%
  ggplot(aes(x = K, y = `Normalized Stat.`, col = method, linetype = method)) +
  geom_line(size=1.1)+
  geom_point(size=2.5)+
  scale_x_continuous(breaks = 1:13) +
  theme_bw() + 
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position="bottom",
    legend.text = ggplot2::element_text(size=13),
    text = ggplot2::element_text(size=22),
    plot.title = element_text(hjust=0.98, vjust=-6),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) + 
  guides(colour=guide_legend(keywidth = 1.5, keyheight = 0.5)) + 
  ggtitle(net_name)+
  theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"))

ggsave(file.path("FBoutput", paste0(net_name, "_statplot.pdf")), width=6, height=5)
saveRDS(out, file.path("FBouput", paste0(net_name, "_norm_stats.rds")))