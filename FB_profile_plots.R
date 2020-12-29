# generate profile plots for facebook network

library(nett)
library(tidyverse)
library(igraph)
library(doParallel)
library(KRLS2)
source("plot_smooth_profile.R")


if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))

college_names <- names(Alist)

# -------------------------------------------------------------------------
# networks in the main text
net_id = which(college_names == "Bucknell39")
net_id = which(college_names == "Oklahoma97")
net_id = which(college_names == "Northeastern19")
net_id = which(college_names == "Maryland58")
net_id = which(college_names == "UCF52")
net_id = which(college_names == "UC61")
net_id = which(college_names == "Stanford3")
net_id = which(college_names == "Harvard1")
net_id = which(college_names == "Wellesley22")

# networks in the appendix 

net_id = which(college_names == "Vassar85")
net_id = which(college_names == "Lehigh96")
net_id = which(college_names == "Cal65")
net_id = which(college_names == "GWU54")
net_id = which(college_names == "UCSB37")
net_id = which(college_names == "Temple83")
net_id = which(college_names == "Baylor93")
net_id = which(college_names == "Swarthmore42")
net_id = which(college_names == "Santa74")
net_id = which(college_names == "Notre Dame57")
net_id = which(college_names == "UCLA26")
net_id = which(college_names == "FSU53")
net_id = which(college_names == "Santa74")
net_id = which(college_names == "Texas80")
net_id = which(college_names == "Michigan23")
net_id = which(college_names == "Columbia2")
net_id = which(college_names == "Yale4")
net_id = which(college_names == "Caltech36")


result_folder = "FBoutput"
nrep = 500

net_name = college_names[net_id]
tag = paste0(net_name, "_SNAC+_smooth_",nrep)
dat_fname = file.path(result_folder, paste0(tag, ".rda"))

if (file.exists(dat_fname)) {
  load(dat_fname)
} else {
  g = graph_from_adjacency_matrix(Alist[[net_id]], mode = "undirected")
  g = extract_largest_cc(g)
  A = extract_low_deg_comp(g, verb=T) %>% as_adj()
  tstat = snac_resample(A, nrep = nrep, ncores = 16)
  save(tstat, file = dat_fname)
}
# plot_smooth_profile(tstat, net_name, krr=T, trunc_type = "none", spar=0.3, plot_null_spar = T)
plot_smooth_profile(tstat, net_name, krr=F, trunc_type = "none", spar=0.3, plot_null_spar = T)
ggsave(file.path(result_folder, paste0(tag, ".pdf")), width=6, height=5)

