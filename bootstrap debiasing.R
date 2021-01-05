# bootstrap debiasing
boot_nac <- function(A, z, K, nrep = 10, dcsbm_boot = T) {
  # printf('nboot = %d ', nrep)
  # bootstrap with SBM
  out <- nett::estim_dcsbm(A, z)
  B <- out$B
  # bootstrap with DCSBM
  # theta <- out$theta
  boot_samp <-sapply(1:nrep, function(i) {
    Aboot <- nett::fast_sbm(z, B)
    # Aboot <- nett::sample_dcsbm(z, B, theta = theta)
    zboot <- spec_clust(Aboot, K)
    c(
      nac_test(Aboot, K, z = zboot)$stat,
      snac_test(Aboot, K, z=zboot)$stat,
      adj_spec_test(Aboot, K, zboot)
    )
  })
  means = apply(boot_samp, 1, mean)
  sds = apply(boot_samp, 1, sd)
  list(mean_nac = means[1], mean_snac = means[2], mean_spec = means[3],
       sd_nac = sds[1], sd_snac = sds[2], sd_spec = sds[3])
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

model_select <- function(A, labels, alpha=1e-4, bootrep=5, Kmin = 1, Kmax) {
  Kmax = ncol(labels)-1
  n = nrow(A)
  Tnac  = sapply(Kmin:Kmax, function(k) nac_test(A, k, z = labels[,k-Kmin+1], y = labels[,k-Kmin+2])$stat )
  Tsnac = sapply(Kmin:Kmax, function(k) snac_test(A, k, z = labels[,k-Kmin+1])$stat )
  Tspec = sapply(Kmin:Kmax, function(k) adj_spec_test(A, k, labels[,k-Kmin+1]))
  bic = sapply(Kmin:Kmax, function(k) eval_dcsbm_bic(A, z = labels[ , k-Kmin+1], k, poi = T))
  
  boots <- do.call(rbind, lapply(Kmin:Kmax, function(k) data.frame(boot_nac(A, z = labels[,k-Kmin+1], K = k, nrep=bootrep))))
  Tnac_de = standardize(Tnac, boots$mean_nac, boots$sd_nac)
  Tsnac_de = standardize(Tsnac, boots$mean_snac, boots$sd_snac)
  Tspec_de = standardize(Tspec, boots$mean_spec, boots$sd_spec)
  
  th_norm = qnorm(alpha, lower.tail = F)
  th_tw <- RMTstat::qtw(alpha, lower.tail = F)
  stats <- list(snac_stat = Tsnac, 
                snac_boot_stat = Tsnac_de,
                nac_stat = Tnac,
                nac_boot_stat = Tnac_de,
                as_stat = Tspec,
                as_boot_stat = Tspec_de
                )
  result <- list(snac = decide_K(Tsnac, th_norm)+Kmin-1, 
                 snac_boot = decide_K(Tsnac_de, th_norm)+Kmin-1, 
                 nac_boot = decide_K(Tnac_de, th_norm)+Kmin-1,
                 as_boot = decide_K(Tspec_de, th_tw)+Kmin-1,
                 LRBIC = which.max(bic)+Kmin-1
                 )
  return(c(result, stat))
}