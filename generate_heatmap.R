library(nett)
library(pheatmap)
n <- 300
Ktru <- 3 
z = c(rep(1,n/3), rep(2, n/3), rep(3, n/3))

beta <- 0.1
Pi <- rep(1, Ktru)/Ktru # distribution of cluster labels
lambda <- 10
B <- pp_conn(n, beta, lambda, Pi)$B
A <- fast_sbm(z, B)

image(A)
Z <- label_vec2mat(z, Ktru)
P <- Z %*% B %*% t(Z)

L = 2
L = 3
y = spec_clust(A, L)
compute_confusion_matrix(z, y)
# y = sample(1:L, size = n, replace = T)
Y = label_vec2mat(y, L)
X = P %*% Y

# generate heatmap
library(RColorBrewer)
col1 <- brewer.pal(12, "Set3")
mydf = data.frame(row.names =  1:n, 
                  cluster = as.character(z))
mymat = X/rowSums(X)
colnames(mymat) = 1:L
rownames(mymat) = 1:n
pheatmap(mymat, cluster_cols = F, cluster_rows = F,
         gaps_row = c(100, 200), border_color = "black",
         annotation_row = mydf, show_rownames = F)


