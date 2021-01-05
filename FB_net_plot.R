# Network plot 
if (!exists("Alist")) Alist <- readRDS(file.path("data","Alist.rds"))
college_names <- names(Alist)
net_id = which(college_names == "Caltech36")
# net_id = which(college_names == "Northeastern19")
# net_id = which(college_names == "Harvard1")
# net_id = which(college_names == "Stanford3")
g = graph_from_adjacency_matrix(Alist[[net_id]], mode = "undirected")
g = extract_largest_cc(g)
net_name = college_names[net_id]

degs = degree(g)
deg_prec = 0.75
g2 = g %>% 
  induced_subgraph(degs <= quantile(degs, deg_prec))%>%
  extract_largest_cc()

# leave out nodes with low degrees
# havard > 10
# northeastern >2
# stanford > 3
# maryland > 0
degs2 <- degree(g2)
idx <- degs2 > 0
g3 <- g2 %>% 
  induced_subgraph(idx)%>%
  extract_largest_cc()

K = 3
dev.off()
A <- as_adjacency_matrix(g3)
zh = spec_clust(A, K, tau = 0.25)
zh2 = nett::fast_cpl(A, K, ilabels = zh, niter = 20)
compute_mutual_info(zh2,zh)

net_tag = sprintf("%s_%dpct", net_name, round(100*deg_prec))
coord = NULL
coord_filename = file.path("data", sprintf("%s_coord.rds", net_tag))
if (file.exists(coord_filename)) {
  coord = readRDS(coord_filename)
} 
par(mar = c(0,0,0,0))

out = plot_net(g3, community = zh2, vertex_alpha=.55, coord = coord)
# title(net_name, cex.main = 3)
# text(0.1,1, sprintf("(degrees < %d%%)", round(100*deg_prec)), cex=2)

table(zh2)
dev.copy(pdf, file.path("FBoutput", sprintf("%s_K%d_netplot.pdf", net_name, K)))
dev.off()

saveRDS(out$coord, coord_filename)
