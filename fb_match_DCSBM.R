library(fastRG) # used to generate Poisson DCSBM
library(dplyr)
library(nett)
if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))
Ktru <- 3

deg_match_res <- do.call(bind_rows, mclapply(1:100, function(net_id) {
  A_fb <- Alist[[net_id]]
  deg <- rowSums(A_fb)
  ave_deg <- mean(deg)
  n <- dim(A_fb)[1]
  
  # use degree as theta 
  theta1 <- deg
  B1 = pp_conn(n, oir = 0.1, lambda = ave_deg,
               pri = rep(1,Ktru), theta = theta1)$B
  # Bernoulli generation
  z1_ber <- sample(1:Ktru, n, replace = T)
  A1_ber <- sample_dcsbm(z1_ber, B1, theta1)
  deg1_ber <- rowSums(A1_ber)
  # Poisson generation
  model1 <- fastRG::dcsbm(theta = theta1, B = B1, expected_degree = ave_deg)
  levels(model1$z) <- 1:Ktru
  z1_poi <- model1$z
  A1_poi <- fastRG::sample_sparse(model1)
  deg1_poi <- rowSums(A1_poi)
  
  # estimate pareto distribution parameters
  deg_trunc = deg[deg >= 2] # do truncation to avoid large variance on theta
  # deg_trunc = deg
  shape_h = min(deg_trunc)
  location_h = length(deg_trunc)/(sum(log(deg_trunc)) - length(deg_trunc)*log(shape_h))
  theta2 <- rpareto(n, location_h, shape_h)
  B2 = pp_conn(n, oir = 0.1, lambda = ave_deg,
               pri = rep(1,Ktru), theta = theta2)$B
  #Bernoulli generation
  z2_ber <- z1_ber
  A2_ber <- nett::sample_dcsbm(z2_ber, B2, theta2) 
  deg2_ber <- rowSums(A2_ber)
  # Poisson generation
  model2 <- fastRG::dcsbm(theta = theta2, B = B2, expected_degree = ave_deg)
  levels(model2$z) <- 1:Ktru
  z2_poi <- model2$z
  A2_poi <- fastRG::sample_sparse(model2)
  deg2_poi <- rowSums(A2_poi)
  
  # results 
  tibble::tribble(
    ~method, ~KS_pvalue, ~ave_deg_diff, ~SNAC_zh, ~SNAC_z, 
    "direct_ber", ks.test(deg, deg1_ber)$p.value, abs(mean(deg1_ber) - ave_deg), 
    snac_test(A1_ber, Ktru)$stat, snac_test(A1_ber, Ktru, z = z1_ber)$stat,
    
    "direct_poi", ks.test(deg, deg1_poi)$p.value, abs(mean(deg1_poi) - ave_deg),
    snac_test(A1_poi, Ktru)$stat, snac_test(A1_poi, Ktru, z = z1_poi)$stat,
    
    "est_ber", ks.test(deg, deg2_ber)$p.value, abs(mean(deg2_ber) - ave_deg),
    snac_test(A2_ber, Ktru)$stat, snac_test(A2_ber, Ktru, z = z2_ber)$stat,
    
    "est_poi", ks.test(deg, deg2_poi)$p.value, abs(mean(deg2_poi) - ave_deg),
    snac_test(A2_poi, Ktru)$stat, snac_test(A2_poi, Ktru, z = z2_poi)$stat,
  )
}, mc.cores = 10))

# generate empirical histograms with zh
p1 = deg_match_res %>%
  filter(method == "direct_ber")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p2 = deg_match_res %>%
  filter(method == "direct_poi")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p3 = deg_match_res %>%
  filter(method == "est_ber")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p4 = deg_match_res %>%
  filter(method == "est_poi")%>%
  ggplot(aes(x = SNAC_zh)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

ggpubr::ggarrange(p1, p2, p3, p4, 
                  labels = c("direct_ber", "direct_poi", "est_ber", "est_poi"),
                  ncol = 2, nrow = 2)

# generate empirical histograms with z
p1z = deg_match_res %>%
  filter(method == "direct_ber")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p2z = deg_match_res %>%
  filter(method == "direct_poi")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p3z = deg_match_res %>%
  filter(method == "est_ber")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

p4z = deg_match_res %>%
  filter(method == "est_poi")%>%
  ggplot(aes(x = SNAC_z)) + 
  theme_bw() +  xlab("SNAC")+
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "white", bins = 15)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 2)+
  theme(
    text = element_text(size=26)
  )

ggpubr::ggarrange(p1z, p2z, p3z, p4z, 
                  labels = c("direct_ber", "direct_poi", "est_ber", "est_poi"),
                  ncol = 2, nrow = 2)

# KS test to see the similarity of degree distribution 
# of FB network and DCSBM
nrow(deg_match_res %>% filter(method == "direct_ber" & 
                                KS_pvalue >= 0.05))/100
nrow(deg_match_res %>% filter(method == "direct_poi" & 
                                KS_pvalue >= 0.05))/100
# using FB degree as theta can get DCSBM's degree close to FB's  
nrow(deg_match_res %>% filter(method == "est_ber" & 
                                KS_pvalue >= 0.05))/100
nrow(deg_match_res %>% filter(method == "est_poi" & 
                                KS_pvalue >= 0.05))/100

deg_hmean <- sapply(Alist, function(A) 1/mean(1/rowSums(A)))
summary(deg_hmean)
