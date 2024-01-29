### Heatmap spotting an outlier

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading distance matrix as an input file obtained from context specific GEMs created

data1= read.delim('/home/shweta/Desktop/ProjectChildren/GEM/DistanceMatrix - Sheet1.tsv',header = T,row.names = 1)
col_fun = colorRamp2(c(0.8, 0.97, 0.98, 0.985, 0.99, 1), c("white", "white", "#ff0000","#cc0000", "#990000", "#660000"))

# Labelling the outlier
ha = rowAnnotation(foo = anno_mark(at = c(32), 
                                   labels = ("CA104SB") ))

#saving the heatmap to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/Heatmap1.pdf", width = 7, height = 7)
Heatmap(as.matrix(data1),right_annotation = ha, row_names_gp = gpar(fontsize = 4), row_dend_width = unit(3, "cm"), column_dend_height = unit(4, "cm") ,show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(9, "cm"),width  = unit(9, "cm"), name = "Hamming Distance", heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"), title_gp = gpar(fontsize = 8), at = c(0.8, 0.97, 0.98, 0.99, 1), labels_gp = gpar(fontsize = 8)))

dev.off()


######################################################

### PCA plot showing an outlier

# Importing the required libraries

library(PCAtools)
library(ggrepel)

# Loading Reaction Binary (1 and 0 as presence or absence of rxns) file as an input file obtained from context specific GEMs created and Basic meta data for performing PCA

count=read.delim("/home/shweta/Desktop/ProjectChildren/GEM/ReactionBinaryTable - Sheet1.tsv",row.names = 1,check.names = FALSE)
meta=read.delim("/home/shweta/Desktop/ProjectChildren/BasicMetadata.txt", header=TRUE, row.names = 1, check.names = FALSE)

p <- pca(count, metadata = meta, removeVar = 0.1)


# str(p)
#p$loadings[1:5,1:5]

# Writing the data in rotated attribute of PCA results to a File

data1=write.table(p$rotated,file="/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated.txt",sep="\t",col.names = NA,quote = FALSE)

# Reading the file in which an additional column showing label was added
pca_data=read.delim("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated1.txt",row.names = 1,check.names = FALSE)

# Plotting a PCA plot and saving in a pdf 
pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA.pdf")


ggplot(pca_data, aes(x=PC1, y=PC2)) + geom_point(size=3, fill="#FFC0CB", colour="#5d3954", shape=21)+
  labs(x="PC1, 10.38% variance",y="PC2, 5.65% variance")+theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
                                                               legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  geom_text_repel(aes(label = pca_data$Label), box.padding = 0.5)+ guides(shape = guide_legend(nrow = 1, override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


dev.off()



######################################################


### PCA plot without outlier


# Importing the required libraries

library(PCAtools)
library(ggrepel)

# Loading Reaction Binary (1 and 0 as presence or absence of rxns) file as an input file obtained from context specific GEMs created and Basic meta data for performing PCA

count=read.delim("/home/shweta/Desktop/ProjectChildren/GEM/ReactionBinaryTable_without_outlier.tsv",row.names = 1,check.names = FALSE)
meta=read.delim("/home/shweta/Desktop/ProjectChildren/BasicMetadata_without_outlier.txt", header=TRUE, row.names = 1, check.names = FALSE)

p <- pca(count, metadata = meta, removeVar = 0.1)


# str(p)
#p$loadings[1:5,1:5]

# Writing the data in rotated attribute of PCA results to a File

data1=write.table(p$rotated,file="/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated_without_outlier.txt",sep="\t",col.names = NA,quote = FALSE)

# Reading the file
pca_data=read.delim("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated_without_outlier.txt",row.names = 1,check.names = FALSE)

# Plotting a PCA plot and saving in a pdf 
pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_without_outlier.pdf")


ggplot(pca_data, aes(x=PC1, y=PC2)) + geom_point(size=3, fill="#bae5ee", colour="#0e5464", shape=21)+
  labs(x="PC1, 10.38% variance",y="PC2, 5.65% variance")+theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
                                                               legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1, override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))


dev.off()


######################################################

### K-means with 4 clusters

# Importing the required libraries

library(tidyverse)
library(cluster)
library(factoextra)
library(openxlsx)

# Loading Reaction Binary (1 and 0 as presence or absence of rxns) file as an input file obtained from context specific GEMs created

data3=read.delim("/home/shweta/Desktop/ProjectChildren/GEM/ReactionBinaryTable_without_outlier.tsv",header=TRUE,row.names = 1)

# Creating transpose of the matrix and writing it to the file
data3_trans=t(data3)

write.table((data3_trans),file="/home/shweta/Desktop/ProjectChildren/GEM/Transpose_binary_RxnSample_without_outlier.txt", quote = FALSE,col.names = NA)

# Kmeans clustering using 4 as the given number of cluster

set.seed(123)
final <- kmeans(data3_trans, 4, nstart = 25)

# Storing the Kmeans results
final$cluster

kmeans_cluster=as.data.frame(final$cluster)

write.table((kmeans_cluster),file="/home/shweta/Desktop/ProjectChildren/R_analysis/kmeans_with_4_cluster_without_outlier.txt", quote = FALSE,col.names = NA)



####

# Generating a PCA plot showing 4 kmeans clusters and storing in a pdf

data = read.delim("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated_without_outlier.txt", row.names = 1)

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/Kmeans_4clusters_without_outlier.pdf")
ggplot(data, aes(PC1,PC2,color=Kmeans_clusters))+geom_point(size=3,aes(fill=Kmeans_clusters), shape=21)+ geom_density_2d() +
  scale_color_manual(values=c("Cluster_1"="#7f4f34","Cluster_2"="#93866c", "Cluster_3"="#603642", "Cluster_4"="#65621a"))+
  scale_fill_manual(values=c("Cluster_1"="#ff9f68","Cluster_2"="#f6e0b5", "Cluster_3"="#c06c84", "Cluster_4"="#dfdc85"))+
  labs(x="PC1, 10.38% variance",y="PC2, 5.65% variance")+
  theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))

dev.off()



######################################################

### PCA plot showing the color scaling according to Location (Urban and Rural)

# Importing the required libraries
library(ggplot2)

# Loading an input file obtained from PCA results

data = read.delim("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated_without_outlier.txt", row.names = 1)

# Generating a PCA plot showing Location and storing in a pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/Location_by_Kmeans_without_outlier.pdf")
ggplot(data, aes(PC1,PC2,color=Location))+geom_point(size=3,aes(fill=Location), shape=21)+
  geom_density_2d() +
  scale_color_manual(values=c("Rural"="#725600","Urban"="#003500"))+
  scale_fill_manual(values=c("Rural"="#bf9000","Urban"="#005900"))+
  labs(x="PC1, 10.38% variance",y="PC2, 5.65% variance")+
  theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))

dev.off()



######################################################

### PCA plot showing 4 clusters with location

data = read.delim("/home/shweta/Desktop/ProjectChildren/R_analysis/PCA_rotated_without_outlier.txt", row.names = 1)

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/Kmeans_4clusters_and_location_without_outlier.pdf")
ggplot(data, aes(PC1,PC2,color=Kmeans_clusters))+geom_point(size=2.5,aes(fill=Kmeans_clusters, shape=Location))+ geom_density_2d() +
  scale_color_manual(values=c("Cluster_1"="#7f4f34","Cluster_2"="#93866c", "Cluster_3"="#603642", "Cluster_4"="#65621a"))+
  scale_fill_manual(values=c("Cluster_1"="#ff9f68","Cluster_2"="#f6e0b5", "Cluster_3"="#c06c84", "Cluster_4"="#dfdc85"))+
  scale_shape_manual(values = c("Rural"=15,"Urban"=16))+
  labs(x="PC1, 10.38% variance",y="PC2, 5.65% variance")+
  theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))

dev.off()




######################################################

### Heatmap showing 4 Clusters final

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading distance matrix as an input file obtained from context specific GEMs created

data1= read.delim('/home/shweta/Desktop/ProjectChildren/GEM/DistanceMatrix - Sheet1_without_outlier.tsv',header = T,row.names = 1)

# Setting a color scale for distance values

col_fun = colorRamp2(c(0.8, 0.9, 0.96, 0.97, 0.98, 0.99, 1), c("yellow", "white", "#fbd4ac", "#ee902e", "#b73a58", "#fcf3e2", "#005900"))

# Loading the basic metadata file for giving the annotations on heatmap

chrt=read.csv("/home/shweta/Desktop/ProjectChildren/BasicMetadata_without_outlier.txt",row.names = 1)

chrt_anno = HeatmapAnnotation(df = chrt, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                              col=list(Location=c(Rural="#bf9000",Urban="#005900"),
                                       Age=colorRamp2(c(10,20,30,40), c("#bae5ee","#1ca9c9","#1687a0","#0e5464")),
                                       Gender=c("Male"="#5d3954","Female"="#FFC0CB")))

# Generating a heatmap and storing it to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/Heatmap_GEM_final.pdf", width = 10, height = 7)
Heatmap(as.matrix(data1), row_dend_width = unit(1.5, "cm"), column_dend_height = unit(2, "cm") ,show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(10, "cm"),width  = unit(10, "cm"), name = "Similarity", top_annotation=chrt_anno, heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"), title_gp = gpar(fontsize = 8), at = c(0.8, 0.97, 0.98, 0.99, 1), labels_gp = gpar(fontsize = 8)))

dev.off()




######################################################


