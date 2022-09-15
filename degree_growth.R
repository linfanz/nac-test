# degree growth  
library(ggplot2)

if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds")) 

num_nodes <- sapply(Alist, nrow)

idSet <- which(num_nodes > 8500 & num_nodes < 10000) 
num_pts <- floor(min(num_nodes[idSet])/100)-1
res <- data.frame(net_id = factor(), n = integer(), hd = double())

for (net_id in idSet){
  A <- Alist[[net_id]]
  n <- nrow(A)
  nodes_0 <- sample(n, 100)
  nodes_i <- nodes_0
  hd <- c()
  ns <- c()
  A_graph = graph_from_adjacency_matrix(A, mode="undirected")
  A_graph = set.vertex.attribute(A_graph, "name", value=as.character(1:n))
  
  for (i in 1:num_pts){
    A_sub <- induced_subgraph(A_graph, nodes_i)
    degs <- degree(A_sub)
    #extract original set node degrees 
    if (i == 1){
      degs <- degs[as.character(nodes_0)]
      degs <- degs[degs > 0]
      nodes_0 <- names(degs)
    }else{
      degs <- degs[nodes_0]
    }
    ns[i] <- length(nodes_i)
    # degs <- degs[as.character(nodes_0)] 
    # degs <- degs[degs > 0] # leave out 0 degree
    hd[i] <- 1/mean(1/degs) 
    
    ## add more nodes 
    nodes_i <- c(nodes_i, sample((1:n)[-nodes_i], 100))
  }
  res <- rbind(res, data.frame(net_id = net_id, n = ns, hd = hd))
}


ggplot(res, aes(x = n, y = hd, color = names(Alist)[net_id])) +
  theme_bw() + ylab("h(d)")+
  geom_line(size = 2)+
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.23, 0.8),
        legend.text = element_text(size=25),
        text = element_text(size=26),
        legend.key.width = unit(4, "line"))
# ggsave("degree_growth.pdf", width = 8.22, height = 6.94)

