load("../../Data/PPI_mat.RData")
load("../../Data/PPI.RData")
write.csv(PPI, file="PPI_cytoscape.csv", quote=F)

df<-data.frame(genes=unique(hits, driver_genes), driver=ifelse(unique(hits, driver_genes) %in% driver_genes, "Yes", "No"))
write.csv(df, file="Hits_driver_anno.csv")

load("/Volumes/home/VUS/VUS/results/20220208/hits.RData")

PPI_sel<-PPI_mat[which(rownames(PPI_mat) %in% driver_genes), 
                 which(colnames(PPI_mat) %in% setdiff(hits, driver_genes))]


rand_sum<-c()
set.seed(193758)
for(iter in 1:1000){
rand_genes<-sample(setdiff(1:ncol(PPI_mat), which(colnames(PPI_mat) %in% driver_genes)), length(which(colnames(PPI_mat) %in% setdiff(hits, driver_genes))))
PPI_rand<-PPI_mat[which(rownames(PPI_mat) %in% driver_genes), 
                 rand_genes]
rand_sum<-c(sum(PPI_rand), rand_sum)
}

sum(rand_sum>sum(PPI_sel))

df<-data.frame(First_neighbours=rand_sum)

pdf("figures/PPI_neigghbours.pdf",5,4)
ggplot(df, aes(x=First_neighbours))+geom_density()+ geom_vline(xintercept = sum(PPI_sel), linetype="dashed", 
                                                               color = "red", size=1)+theme_classic()

dev.off()
