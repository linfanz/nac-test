# bootstrap debiasing
# simulate from DCSBM
# add quantile method back
boot_all_nac <- function(A, z, K, nrep = 10, alpha = 1e-6) {
  out <- nett::estim_dcsbm(A, z)
  B <- out$B
  theta <- out$theta
  boot_samp <- sapply(1:nrep, function(i) {
    Aboot <- nett::fast_sbm(z, B)
    AbootDC <- nett::sample_dcsbm(z, B, theta = theta)
    ### reestimate labels 
    zboot <- spec_clust(Aboot, K)
    c(
      nac_test(Aboot, K, z = zboot)$stat,
      snac_test(Aboot, K, z = zboot)$stat,
      # DCSBM generation
      nac_test(AbootDC, K, z = zboot)$stat,
      snac_test(AbootDC, K, z = zboot)$stat
    )
    
    ### use true labels
    # c(
    #   nac_test(Aboot, K, z = z)$stat,
    #   snac_test(Aboot, K, z = z)$stat,
    #   # DCSBM generation
    #   nac_test(AbootDC, K, z = z)$stat,
    #   snac_test(AbootDC, K, z = z)$stat
    # )
  })
  means = apply(boot_samp, 1, mean)
  sds = apply(boot_samp, 1, sd)
  ths = apply(boot_samp, 1, quantile, 1-alpha)
  list(mean_nac = means[1], mean_snac = means[2], 
       sd_nac = sds[1], sd_snac = sds[2], 
       th_nac = ths[1], th_snac = ths[2], 
       # DCSBM generation 
       mean_nac_dc = means[3], mean_snac_dc = means[4], 
       sd_nac_dc = sds[3], sd_snac_dc = sds[4], 
       th_nac_dc = ths[3], th_snac_dc = ths[4]
       )
}

standardize <- function(x,mu,sig) (x - mu)/sig
# Function to decide K based on sequential test statistics
decide_K <- function(Tstat_seq, th) {
  inrange <- Tstat_seq < th
  K <- which(inrange)[1]
  # Kmax is not enough, we return Kmax for consistency
  if (is.na(K)) K <- length(Tstat_seq) 
  return(K)
}

nac_model_select <- function(A, labels, alpha=1e-4, bootrep=5, Kmin = 1, Kmax) {
  Kmax = ncol(labels)-1
  n = nrow(A)
  Tnac  = sapply(Kmin:Kmax, function(k) nac_test(A, k, z = labels[,k-Kmin+1], y = labels[,k-Kmin+2])$stat )
  Tsnac = sapply(Kmin:Kmax, function(k) snac_test(A, k, z = labels[,k-Kmin+1])$stat )
  
  boots <- do.call(rbind, lapply(Kmin:Kmax, function(k) 
    data.frame(boot_all_nac(A, z = labels[,k-Kmin+1], K = k, nrep=bootrep, alpha = alpha))))
  Tnac_de = standardize(Tnac, boots$mean_nac, boots$sd_nac)
  Tsnac_de = standardize(Tsnac, boots$mean_snac, boots$sd_snac)
  # DCSBM generation bootstrap stat
  Tnac_de_dc = standardize(Tnac, boots$mean_nac_dc, boots$sd_nac_dc)
  Tsnac_de_dc = standardize(Tsnac, boots$mean_snac_dc, boots$sd_snac_dc)
  
  th_norm = qnorm(alpha, lower.tail = F)
  # stats <- list(snac_stat = Tsnac, 
  #               snac_boot_stat = Tsnac_de,
  #               nac_stat = Tnac,
  #               nac_boot_stat = Tnac_de
  # )
  result <- list(snac = decide_K(Tsnac, th_norm)+Kmin-1, 
                 ## SBM 
                 # debiasing
                 snac_boot_de = decide_K(Tsnac_de, th_norm)+Kmin-1, 
                 nac_boot_de = decide_K(Tnac_de, th_norm)+Kmin-1,
                 # quantile method 
                 snac_boot_qu = decide_K(Tsnac, boots$th_snac)+Kmin-1, 
                 nac_boot_qu = decide_K(Tnac, boots$th_nac)+Kmin-1,
                 
                 ## DCSBM
                 # debiasing
                 snac_dc_boot_de = decide_K(Tsnac_de_dc, th_norm)+Kmin-1, 
                 nac_dc_boot_de = decide_K(Tnac_de_dc, th_norm)+Kmin-1,
                 # quantile method
                 snac_dc_boot_qu = decide_K(Tsnac, boots$th_snac_dc)+Kmin-1, 
                 nac_dc_boot_qu = decide_K(Tnac, boots$th_nac_dc)+Kmin-1
  )
  # return(c(result, stat))
  return(result)
}