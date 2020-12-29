snac_resample <- function(A, nrep = 20, Kmin = 1, Kmax = 13, ncores = 10, seed = 1234) {
  Ks = Kmin:Kmax
  
  cl <- parallel::makeForkCluster(ncores)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed)
  
  labels = sapply(Kmin:(Kmax+1), function(k) nett::spec_clust(A, k))
  Tstat = foreach(t = 1:nrep) %dopar% {
    tibble(
      itr = rep(t, length(Ks)), 
      K = Ks,
      value = sapply(Ks, function(k) snac_test(A, k, labels[ , k-Kmin+1])$stat)
    )
  } %>% bind_rows()
  
  parallel::stopCluster(cl)
  Tstat
}

# fp = f', fpp = f''
detect_dip = function(fp){
  which(fp > 0)[1]
}
compute_curvature = function(fp, fpp) {
  abs(fpp) / ((1 + fp^2)^(1.5))
}

detect_elbow = function(curvature) {
  which.max(curvature)[1]
}

fit_ss = function(x, y, xx, spar = NULL, trunc_type = "none") {
  trunc_fun = switch(trunc_type, "none" = function(x) signif(x, 2), 
                     "floor" = floor, "ceil" = ceiling, "round" = round)
  
  fitted_ss = smooth.spline(x, y, spar = spar)
  fxx_ss = predict(fitted_ss, xx, deriv = 0)$y
  # first derivative
  fp_ss = predict(fitted_ss, xx, deriv = 1)$y
  # second derivative
  fpp_ss = predict(fitted_ss, xx, deriv = 2)$y
  # curvature 
  # curvature_ss = compute_curvature(fp_ss, fpp_ss)
  
  # plot(xx, fp_ss)
  # plot(xx, fpp_ss)
  # plot(xx, fxx_ss)
  # plot(xx, curvature_ss, log = "y")
  # plot(xx, curvature_ss)
  
  list(dip = trunc_fun( xx[detect_dip(fp_ss)]), 
       # max second derivative
       elbow1 = trunc_fun( xx[detect_elbow(fpp_ss)] ),
       # elbow2 = trunc_fun( xx[detect_elbow(curvature_ss)]),
       fxx = fxx_ss)
}

plot_smooth_profile = function(
  tstat, net_name, krr=T, trunc_type = "none", 
  spar=NULL, plot_null_spar = F,
  alpha = 0.3 # transparency
) {
  trunc_fun = switch(trunc_type, "none" = function(x) signif(x, 2), 
                     "floor" = floor, "ceil" = ceiling, "round" = round)
  
  y = tstat$value 
  x = tstat$K
  Kmax = max(x)
  xx <- seq(1, Kmax, length.out = 1000)
  
  ss_res1 = fit_ss(x, y, xx, spar = spar, trunc_type = trunc_type)
  
  # fxx_ss <- predict(smooth.spline(x, y), xx)$y
  # # Kh_ss = floor(xx[which(diff(fxx_ss)>0)[1]])
  # Kh_ss = trunc_fun(xx[which(diff(fxx_ss)>0)[1]])
  
  # fx <- predict(smooth.spline(x, y), 1:Kmax)$y
  # (1:Kmax)[which.min(abs(diff(fx)))[1]]
  
  # Kh1 = floor(xx[which.min(abs(diff(fxx_ss)))[1]])
  # Kh2 = floor(xx[which(diff(fxx_ss)>0)[1]])
  # Kh = min(Kh1, Kh2)
  if (krr) {
    # to install this version of KRLS2, run 
    # devtools::install_github("linfanz/KRLS@apprx")
    require(KRLS2)
    fxx_krr <- predict(krls_basic(y, x, b = 3, folds = 10), xx)
    # dx = xx[2] - xx[1]
    dx = 1
    fp_krr = diff(fxx_krr)/dx
    fpp_krr = diff(fp_krr)/dx
    curvature_krr = compute_curvature(fp_krr[-length(fp_krr)], fpp_krr)
    # Kh_krr = trunc_fun(xx[which(diff(fxx_krr)>0)[1]])
    dip_krr = trunc_fun(xx[detect_dip(fp_krr) + 1])  # + 1 since we are using diff
    elbow_krr = trunc_fun(xx[detect_elbow(curvature_krr)+1])
  }
  
  color2 = "#FB674B"
  
  elbow <- ifelse(is.na(ss_res1$elbow1),"NA",ss_res1$elbow1)
  dip <- ifelse(is.na(ss_res1$dip),"NA",ss_res1$dip)
  
  # optK_label = bquote(atop(.(net_name), K %~~% .(ss_res1$elbow1) ~","~.(ss_res1$elbow2)~","~.(ss_res1$dip)))
  # optK_label = c(bquote(.(net_name)), bquote(K %~~% .(ss_res1$elbow1) ~","~.(ss_res1$elbow2)~","~.(ss_res1$dip)))
  # optK_label = c(shQuote(net_name), bquote(K %~~% .(ss_res1$elbow1) ~","~.(ss_res1$elbow2)~","~.(ss_res1$dip)))
  optK_label = c(shQuote(net_name), bquote(K %~~% .(elbow) ~","~.(dip)))
  p = ggplot(tstat, aes(x=K, y=value)) + 
    geom_point(size = 2, color="#6F98C2", alpha = alpha) + 
    geom_line(aes(x=xx, y=fxx), data = tibble(xx=xx, fxx=ss_res1$fxx), size = 1.1, color="black")  +
    scale_x_continuous(breaks = unique(tstat$K)) +
    # ylab("SNAC+") +
    theme_bw(base_size = 20)
  if (krr) {
    color2 = "blue"
    p = p + geom_line(aes(x=xx, y=fxx), data = tibble(xx=xx, fxx=fxx_krr), size = 3, color=color2, linetype = "dashed")
    # annotate("text", x=Kmax-1.5, y=max(tstat$value)*0.88, 
    #          # label = bquote(K %~~% ~.(Kh_krr)), size = 5, color="blue")
    #          label = bquote(K %~~% ~.(elbow_krr))~","~.(dip_krr), size = 5, color="blue")
    optK_label = c(optK_label, bquote(K %~~% .(elbow_krr) ~","~.(dip_krr)))
  }

  if (!is.null(spar) && plot_null_spar & !krr) {
    ss_res2 = fit_ss(x, y, xx, spar = NULL, trunc_type = trunc_type)
    elbow <- ifelse(is.na(ss_res2$elbow1),"NA",ss_res2$elbow1)
    dip <- ifelse(is.na(ss_res2$dip),"NA",ss_res2$dip)
    p = p + 
      # annotate("text", x=Kmax-1.5, y=max(tstat$value)*0.88, 
      #          label = bquote(K %~~% .(ss_res2$elbow1) ~","~.(ss_res2$elbow2)~","~.(ss_res2$dip)), size = 5, color="#FB674B")+
      geom_line(aes(x=xx, y=fxx), data = tibble(xx=xx, fxx=ss_res2$fxx), size = 1.5, color = color2, linetype = "dashed")
    # optK_label =c(optK_label, bquote(K %~~% .(ss_res2$elbow1) ~","~.(ss_res2$elbow2)~","~.(ss_res2$dip)))
    # optK_label =c(optK_label, bquote(K %~~% .(ss_res2$elbow1) ~","~.(ss_res2$dip)))
    optK_label =c(optK_label, bquote(K %~~% .(elbow) ~","~.(dip)))
    
  }
  p + annotate("text", x=Kmax-2.3, y=max(tstat$value)*c(0.98,0.91,0.84),
               # label = bquote(atop(.(net_name), K %~~% ~.(Kh_ss))), size = 5)+
               label = optK_label, size = 7, parse =T, color = c("black", "black", color2))+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "cm")
    )
}





