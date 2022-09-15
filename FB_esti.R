# estimate B using all FB networks
# save the estimated theta!
library(nett)
if (!exists("Alist")) Alist <- readRDS(file.path("FBdata","Alist.rds"))
# if (!exists("Alist_reduced")) Alist_reduced <- readRDS(file.path("FBdata","Alist_75reduced.rds"))

Ktru = 3
block_sums = block_ns = matrix(0, nrow = Ktru, ncol = Ktru)
theta_list = list()
zh_list = list()
for(i in 1:100){
  A = Alist[[i]]
  # A = Alist_reduced[[i]]
  zh = spec_clust(A, Ktru)
  zh_list[[i]] = zh
  block_sums_i = compute_block_sums(A, zh)
  nhs_i = table(zh)
  block_ns_i = nhs_i %*% t(nhs_i) - diag(nhs_i)
  
  block_sums = block_sums +  block_sums_i 
  block_ns = block_ns + block_ns_i
  
  theta_i = estim_dcsbm(A, zh)$theta
  theta_list[[i]] = theta_i
}

B = block_sums/block_ns

saveRDS(theta_list, file = sprintf("FBdata/fb_theta_est_K=%d.rds", Ktru))
saveRDS(B, file = sprintf("FBdata/fb_B_est_K=%d.rds", Ktru))
saveRDS(zh_list, file =  sprintf("FBdata/fb_z_est_K=%d.rds", Ktru))

# saveRDS(theta_list, file = sprintf("FBdata/fb_reduced_theta_est_K=%d.rds", Ktru))
# saveRDS(B, file = sprintf("FBdata/fb_reduced_B_est_K=%d.rds", Ktru))
# saveRDS(zh_list, file =  sprintf("FBdata/fb_reduced_z_est_K=%d.rds", Ktru))



