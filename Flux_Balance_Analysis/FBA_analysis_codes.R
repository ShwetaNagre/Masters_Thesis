### Splitting master data into clusters

# Mixed_2 is a master file obtained after performing FBA. We got 150 rxns id's files and 150 flux values files,
# and the we merged and mixed them to get this master file containing Flux values for all the 150 samples.
# The clusters obtained from Kmeans algorithm was added to the file

uni_rxn_data = read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Mixed_2.txt", row.names = 1, header = TRUE)

# NA values were replaced with 0
uni_rxn_data[is.na(uni_rxn_data)] = 0

# Transposed the matrix
transposed_uni_rxn = t(uni_rxn_data)

transposed_uni_rxn <- as.data.frame(transposed_uni_rxn)

# The samples in the matrix were separated based on the cluster
cluster_data = split(transposed_uni_rxn, transposed_uni_rxn$Clusters)
View(cluster_data)

# All the samples and there fluxes belonging to a particular cluster were stored in the files
C1=t(cluster_data$Cluster_1)

C2=t(cluster_data$Cluster_2)

C3=t(cluster_data$Cluster_3)

C4=t(cluster_data$Cluster_4)

write.table(C1, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_rxns.txt", sep="\t",col.names = NA,quote = FALSE)

write.table(C2, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_rxns.txt", sep="\t",col.names = NA,quote = FALSE)

write.table(C3, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_rxns.txt", sep="\t",col.names = NA,quote = FALSE)

write.table(C4, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_rxns.txt", sep="\t",col.names = NA,quote = FALSE)

#######################################################


### Unique reactions in cluster1

# The C1 file, that contains all the rxns and their fluxes belonging to cluster 1 was used as an input
Cluster1_file=read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_rxns.txt", header = TRUE)

# Reactions with all positive values
all_positive <- Cluster1_file[rowSums(Cluster1_file[, -1] > 0) == ncol(Cluster1_file) -1, ]

# Reactions with all negative values
all_negative <- Cluster1_file[rowSums(Cluster1_file[, -1] < 0) == ncol(Cluster1_file) -1, ]

# Reactions with all zero values
all_zero <- Cluster1_file[rowSums(Cluster1_file[, -1] == 0) == ncol(Cluster1_file) -1, ]

#save the results
C1_positive=write.table(all_positive, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_positive_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C1_negative=write.table(all_negative, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_negative_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C1_zero=write.table(all_zero, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_zero_rxns.txt",sep="\t",col.names = NA,quote = FALSE)


#######################################################


### Unique reactions in cluster2

# The C2 file, that contains all the rxns and their fluxes belonging to cluster 2 was used as an input
Cluster2_file=read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_rxns.txt", header = TRUE)

# Reactions with all positive values
all_positive <- Cluster2_file[rowSums(Cluster2_file[, -1] > 0) == ncol(Cluster2_file) -1, ]

# Reactions with all negative values
all_negative <- Cluster2_file[rowSums(Cluster2_file[, -1] < 0) == ncol(Cluster2_file) -1, ]

# Reactions with all zero values
all_zero <- Cluster2_file[rowSums(Cluster2_file[, -1] == 0) == ncol(Cluster2_file) -1, ]

#save the results
C2_positive=write.table(all_positive, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_positive_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C2_negative=write.table(all_negative, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_negative_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C2_zero=write.table(all_zero, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_zero_rxns.txt",sep="\t",col.names = NA,quote = FALSE)


#######################################################

### Unique reactions in cluster3

# The C3 file, that contains all the rxns and their fluxes belonging to cluster 3 was used as an input
Cluster3_file=read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_rxns.txt", header = TRUE)

# Reactions with all positive values
all_positive <- Cluster3_file[rowSums(Cluster3_file[, -1] > 0) == ncol(Cluster3_file) -1, ]

# Reactions with all negative values
all_negative <- Cluster3_file[rowSums(Cluster3_file[, -1] < 0) == ncol(Cluster3_file) -1, ]

# Reactions with all zero values
all_zero <- Cluster3_file[rowSums(Cluster3_file[, -1] == 0) == ncol(Cluster3_file) -1, ]

#save the results
C3_positive=write.table(all_positive, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_positive_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C3_negative=write.table(all_negative, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_negative_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C3_zero=write.table(all_zero, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_zero_rxns.txt",sep="\t",col.names = NA,quote = FALSE)


#######################################################

### Unique reactions in cluster4

# The C4 file, that contains all the rxns and their fluxes belonging to cluster 4 was used as an input
Cluster4_file=read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_rxns.txt", header = TRUE)

# Reactions with all positive values
all_positive <- Cluster4_file[rowSums(Cluster4_file[, -1] > 0) == ncol(Cluster4_file) -1, ]

# Reactions with all negative values
all_negative <- Cluster4_file[rowSums(Cluster4_file[, -1] < 0) == ncol(Cluster4_file) -1, ]

# Reactions with all zero values
all_zero <- Cluster4_file[rowSums(Cluster4_file[, -1] == 0) == ncol(Cluster4_file) -1, ]

#save the results
C4_positive=write.table(all_positive, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_positive_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C4_negative=write.table(all_negative, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_negative_rxns.txt",sep="\t",col.names = NA,quote = FALSE)

C4_zero=write.table(all_zero, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_zero_rxns.txt",sep="\t",col.names = NA,quote = FALSE)


#######################################################


### C1 merging fluxes

# The C1_all_unique_rxns input file is obtained from Interactive Venn online tool.
# It has the reactions ID's having all positive, all negavtive, and all zero fluxes for the samples
# This file is merged with the master file to fetch the fluxes of all these unique reactions

small_list <- read.table("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_all_unique_rxns.txt", header = TRUE)
big_list <- read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Mixed_2.txt", header = TRUE)

merged_data <- merge(small_list, big_list, by = 'ID')

fetched_reactions <- merged_data$ID
fetched_values <- merged_data[, 2:ncol(merged_data)]

C1_unique <- cbind(fetched_reactions, fetched_values)

# Saving the fetched fluxes in a file
write.table(C1_unique, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_unique_rxns_final.txt",sep="\t",col.names = NA,quote = FALSE)


#############################################################

### C2 merging fluxex

# The C2_all_unique_rxns input file is obtained from Interactive Venn online tool.
# It has the reactions ID's having all positive, all negavtive, and all zero fluxes for the samples
# This file is merged with the master file to fetch the fluxes of all these unique reactions

small_list <- read.table("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_all_unique_rxns.txt", header = TRUE)
big_list <- read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Mixed_2.txt", header = TRUE)

merged_data <- merge(small_list, big_list, by = 'ID')

fetched_reactions <- merged_data$ID
fetched_values <- merged_data[, 2:ncol(merged_data)]

C2_unique <- cbind(fetched_reactions, fetched_values)

# Saving the fetched fluxes in a file
write.table(C2_unique, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_unique_rxns_final.txt",sep="\t",col.names = NA,quote = FALSE)


#############################################################

### C3 merging fluxex

# The C3_all_unique_rxns input file is obtained from Interactive Venn online tool.
# It has the reactions ID's having all positive, all negavtive, and all zero fluxes for the samples
# This file is merged with the master file to fetch the fluxes of all these unique reactions

small_list <- read.table("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_all_unique_rxns.txt", header = TRUE)
big_list <- read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Mixed_2.txt", header = TRUE)

merged_data <- merge(small_list, big_list, by = 'ID')

fetched_reactions <- merged_data$ID
fetched_values <- merged_data[, 2:ncol(merged_data)]

C3_unique <- cbind(fetched_reactions, fetched_values)

# Saving the fetched fluxes in a file
write.table(C3_unique, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_unique_rxns_final.txt",sep="\t",col.names = NA,quote = FALSE)


#############################################################

### C4 merging fluxex

# The C4_all_unique_rxns input file is obtained from Interactive Venn online tool.
# It has the reactions ID's having all positive, all negavtive, and all zero fluxes for the samples
# This file is merged with the master file to fetch the fluxes of all these unique reactions

small_list <- read.table("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_all_unique_rxns.txt", header = TRUE)
big_list <- read.csv("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Mixed_2.txt", header = TRUE)

merged_data <- merge(small_list, big_list, by = 'ID')

fetched_reactions <- merged_data$ID
fetched_values <- merged_data[, 2:ncol(merged_data)]

C4_unique <- cbind(fetched_reactions, fetched_values)

# Saving the fetched fluxes in a file
write.table(C4_unique, file = "/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_unique_rxns_final.txt",sep="\t",col.names = NA,quote = FALSE)




#############################################################

### subsystem count for Cluster1

# The subsystems to which the unique rxns belong were obtained and then counted 
C1_subsys_count = read.delim("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_rxns_subsys_final.txt", row.names = 1)

C1_sub=table(C1_subsys_count$SubSystem)

subsys=write.table(C1_sub,file="/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_subsys_count.txt",sep="\t",col.names = NA,quote = FALSE)


#############################################################

### subsystem count for Cluster2

# The subsystems to which the unique rxns belong were obtained and then counted 
C2_subsys_count = read.delim("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_rxns_subsys_final.txt", row.names = 1)

C2_sub=table(C2_subsys_count$SubSystem)

subsys=write.table(C2_sub,file="/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_subsys_count.txt",sep="\t",col.names = NA,quote = FALSE)



#############################################################

### subsystem count for Cluster3

# The subsystems to which the unique rxns belong were obtained and then counted 
C3_subsys_count = read.delim("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_rxns_subsys_final.txt", row.names = 1)

C3_sub=table(C3_subsys_count$SubSystem)

subsys=write.table(C3_sub,file="/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_subsys_count.txt",sep="\t",col.names = NA,quote = FALSE)



#############################################################

### subsystem count for Cluster4

# The subsystems to which the unique rxns belong were obtained and then counted 
C4_subsys_count = read.delim("/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_rxns_subsys_final.txt", row.names = 1)

C4_sub=table(C4_subsys_count$SubSystem)

subsys=write.table(C4_sub,file="/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_subsys_count.txt",sep="\t",col.names = NA,quote = FALSE)


#############################################################


### Heatmap for Cluster 1 Unique Rxns

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading the input file containing unique reactions and there fluxes for Cluster 1 samples

Uni= read.csv('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_unique_rxns_final.txt',header = T, row.names = 1)
col_fun = colorRamp2(c(-1000, -10, 0, 10, 1000), c("#ffa500", "#FFFF00","white", "#800080","#0000FF"))

Uni[is.na(Uni)] <- 0  # Replace NA with 0

# Basic metadata having the significant parameters is used for the top annotation

chrt=read.csv("/home/shweta/Desktop/ProjectChildren/Basic_metadata_with_signi_parameters.txt", header = T, row.names = 1)

chrt_anno = HeatmapAnnotation(df = chrt, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                              col=list(Location=c(Rural="#bf9000",Urban="#005900"),
                                       Winter_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 9), c("#bae5ee","#1ca9c9","#1687a0","#0e5464")),
                                       Summer_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 12), c("#fcf5b2","#fbee7f","#f9e432","#f8de00")),
                                       Immunisation=c("Complete"="#727472","Incomplete"="#c0d6e4"),
                                       Peanut_Exposure=c("No"="#6c2a28","Yes"="#b59493"),
                                       Regular_Contact_with_Farm_Animals=c("No"="#5d3954","Yes"="#FFC0CB"),
                                       Clusters=c("Cluster_1"="#a23443","Cluster_2"="#c0c0c0", "Cluster_3"="#b47fb4", "Cluster_4"="#cade71")))


# The Subsystems for the unique reactions were labelled as the left annotation

subsys=read.delim('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C1_rxns_subsys_final.txt', row.names = 1, header = T)
subsys=subsys[4:4]


subsys_anno = rowAnnotation(df = subsys, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                            col=list(SubSystem=c("Carnitine shuttle (mitochondrial)" = "#7e6c91",
                                                 "Fatty acid oxidation" = "#00d6be",
                                                 "Omega-3 fatty acid metabolism" = "#ffb999",
                                                 "Others" = "#ab1222",
                                                 "Purine metabolism" = "#cade71",
                                                 "Transport reactions" = "#d8cea1")))


# Generating a heatmap and storing it to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/FBA_analysis/C1_heatmap.pdf", width = 12, height = 7)
Heatmap(as.matrix(Uni), column_dend_height = unit(2, "cm"), row_names_centered = TRUE,
        show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(10, "cm"), width  = unit(10, "cm"),
        name = "Flux_value", top_annotation = chrt_anno, left_annotation = subsys_anno, column_split = chrt$Clusters, row_split = subsys$SubSystem,
        heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                    title_gp = gpar(fontsize = 8), at = c(-1000, -10, 0, 10, 1000), labels_gp = gpar(fontsize = 8)), row_title = " ")

dev.off()





#############################################################


### Heatmap for Cluster 2 Unique Rxns

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading the input file containing unique reactions and there fluxes for Cluster 2 samples

Uni= read.csv('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_unique_rxns_final.txt',header = T, row.names = 1)
col_fun = colorRamp2(c(-1000, -10, 0, 10, 1000), c("#ffa500", "#FFFF00","white", "#800080","#0000FF"))

Uni[is.na(Uni)] <- 0  # Replace NA with 0

# Basic metadata having the significant parameters is used for the top annotation

chrt=read.csv("/home/shweta/Desktop/ProjectChildren/Basic_metadata_with_signi_parameters.txt", header = T, row.names = 1)

chrt_anno = HeatmapAnnotation(df = chrt, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                              col=list(Location=c(Rural="#bf9000",Urban="#005900"),
                                       Winter_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 9), c("#bae5ee","#1ca9c9","#1687a0","#0e5464")),
                                       Summer_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 12), c("#fcf5b2","#fbee7f","#f9e432","#f8de00")),
                                       Immunisation=c("Complete"="#727472","Incomplete"="#c0d6e4"),
                                       Peanut_Exposure=c("No"="#6c2a28","Yes"="#b59493"),
                                       Regular_Contact_with_Farm_Animals=c("No"="#5d3954","Yes"="#FFC0CB"),
                                       Clusters=c("Cluster_1"="#a23443","Cluster_2"="#c0c0c0", "Cluster_3"="#b47fb4", "Cluster_4"="#cade71")))


# The Subsystems for the unique reactions were labelled as the left annotation

subsys=read.delim('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C2_rxns_subsys_final.txt', row.names = 1, header = T)
subsys=subsys[4:4]


subsys_anno = rowAnnotation(df = subsys, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                            col=list(SubSystem=c("Amino acid metabolism" = "#1c2331",
                                                 "Beta oxidation of unsaturated fatty acids (n-9) (mitochondrial)" = "#7e6c91",
                                                 "Fatty acid biosynthesis" = "#00d6be",
                                                 "Fatty acid oxidation" = "#00ffa6",
                                                 "Leukotriene metabolism" = "#ffb999",
                                                 "Omega fatty acid metabolism" = "#f5b41a", "Others" = "#ab1222",
                                                 "Transport reactions" = "#d8cea1")))

# Generating a heatmap and storing it to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/FBA_analysis/C2_heatmap.pdf", width = 12, height = 7)
Heatmap(as.matrix(Uni), column_dend_height = unit(2, "cm"), row_names_centered = TRUE,
        show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(10, "cm"), width  = unit(10, "cm"),
        name = "Flux_value", top_annotation = chrt_anno, left_annotation = subsys_anno, column_split = chrt$Clusters, row_split = subsys$SubSystem,
        heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                    title_gp = gpar(fontsize = 8), at = c(-1000, -10, 0, 10, 1000), labels_gp = gpar(fontsize = 8)), row_title = " ")

dev.off()



#############################################################


### Heatmap for Cluster 3 Unique Rxns

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading the input file containing unique reactions and there fluxes for Cluster 3 samples

Uni= read.csv('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_unique_rxns_final.txt',header = T, row.names = 1)
col_fun = colorRamp2(c(-1000, -10, 0, 10, 1000), c("#ffa500", "#FFFF00","white", "#800080","#0000FF"))

Uni[is.na(Uni)] <- 0  # Replace NA with 0

# Basic metadata having the significant parameters is used for the top annotation

chrt=read.csv("/home/shweta/Desktop/ProjectChildren/Basic_metadata_with_signi_parameters.txt", header = T, row.names = 1)

chrt_anno = HeatmapAnnotation(df = chrt, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                              col=list(Location=c(Rural="#bf9000",Urban="#005900"),
                                       Winter_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 9), c("#bae5ee","#1ca9c9","#1687a0","#0e5464")),
                                       Summer_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 12), c("#fcf5b2","#fbee7f","#f9e432","#f8de00")),
                                       Immunisation=c("Complete"="#727472","Incomplete"="#c0d6e4"),
                                       Peanut_Exposure=c("No"="#6c2a28","Yes"="#b59493"),
                                       Regular_Contact_with_Farm_Animals=c("No"="#5d3954","Yes"="#FFC0CB"),
                                       Clusters=c("Cluster_1"="#a23443","Cluster_2"="#c0c0c0", "Cluster_3"="#b47fb4", "Cluster_4"="#cade71")))


# The Subsystems for the unique reactions were labelled as the left annotation

subsys=read.delim('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C3_rxns_subsys_final.txt', row.names = 1, header = T)
subsys=subsys[4:4]


subsys_anno = rowAnnotation(df = subsys, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                            col=list(SubSystem=c("Carnitine shuttle (ER, cytosolic, mitochondrial)" = "#7e6c91",
                                                 "Fatty acid oxidation" = "#00d6be",
                                                 "Others" = "#ab1222",
                                                 "Transport reactions" = "#d8cea1")))

# Generating a heatmap and storing it to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/FBA_analysis/C3_heatmap.pdf", width = 12, height = 7)
Heatmap(as.matrix(Uni), column_dend_height = unit(2, "cm"), row_names_centered = TRUE,
        show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(10, "cm"), width  = unit(10, "cm"),
        name = "Flux_value", top_annotation = chrt_anno, left_annotation = subsys_anno, column_split = chrt$Clusters, row_split = subsys$SubSystem,
        heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                    title_gp = gpar(fontsize = 8), at = c(-1000, -10, 0, 10, 1000), labels_gp = gpar(fontsize = 8)), row_title = " ")

dev.off()




#############################################################


### Heatmap for Cluster 4 Unique Rxns

# Importing the required libraries

library(ComplexHeatmap)
library(circlize)

# Loading the input file containing unique reactions and there fluxes for Cluster 4 samples

Uni= read.csv('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_unique_rxns_final.txt',header = T, row.names = 1)
col_fun = colorRamp2(c(-1000, -10, 0, 10, 1000), c("#ffa500", "#FFFF00","white", "#800080","#0000FF"))

Uni[is.na(Uni)] <- 0  # Replace NA with 0

# Basic metadata having the significant parameters is used for the top annotation

chrt=read.csv("/home/shweta/Desktop/ProjectChildren/Basic_metadata_with_signi_parameters.txt", header = T, row.names = 1)

chrt_anno = HeatmapAnnotation(df = chrt, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                              col=list(Location=c(Rural="#bf9000",Urban="#005900"),
                                       Winter_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 9), c("#bae5ee","#1ca9c9","#1687a0","#0e5464")),
                                       Summer_Sunlight_Exposure=colorRamp2(c(0, 1, 5, 12), c("#fcf5b2","#fbee7f","#f9e432","#f8de00")),
                                       Immunisation=c("Complete"="#727472","Incomplete"="#c0d6e4"),
                                       Peanut_Exposure=c("No"="#6c2a28","Yes"="#b59493"),
                                       Regular_Contact_with_Farm_Animals=c("No"="#5d3954","Yes"="#FFC0CB"),
                                       Clusters=c("Cluster_1"="#a23443","Cluster_2"="#c0c0c0", "Cluster_3"="#b47fb4", "Cluster_4"="#cade71")))


# The Subsystems for the unique reactions were labelled as the left annotation

subsys=read.delim('/home/shweta/Desktop/ProjectChildren/FBA/All_FBAs/Merged/Unique_rxns/C4_rxns_subsys_final.txt', row.names = 1, header = T)
subsys=subsys[4:4]


subsys_anno = rowAnnotation(df = subsys, simple_anno_size = unit(0.25, "cm"),annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 8), annotation_legend_param  = list(grid_width = unit(0.3, "cm"), grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                            col=list(SubSystem=c("Beta oxidation of fatty acids (mitochondrial)" = "#1c2331",
                                                 "Bile acid biosynthesis" = "#7e6c91",
                                                 "Carnitine shuttle (mitochondrial, ER)" = "#00d6be",
                                                 "Fatty acid activation (cytosolic)" = "#00ffa6",
                                                 "Fatty acid biosynthesis" = "#ffb999",
                                                 "Fatty acid elongation (even-chain)" = "#f7d08a",
                                                 "Fatty acid oxidation" = "#bd707b",
                                                 "Others" = "#ab1222",
                                                 "Purine metabolism" = "#cade71",
                                                 "Pyrimidine metabolism" = "#ff9f33",
                                                 "Transport reactions" = "#d8cea1")))

# Generating a heatmap and storing it to pdf

pdf("/home/shweta/Desktop/ProjectChildren/R_analysis/FBA_analysis/C4_heatmap.pdf", width = 12, height = 7)
Heatmap(as.matrix(Uni), column_dend_height = unit(2, "cm"), row_names_centered = TRUE,
        show_row_names = FALSE, show_column_names = FALSE, col = col_fun, height  = unit(10, "cm"), width  = unit(10, "cm"),
        name = "Flux_value", top_annotation = chrt_anno, left_annotation = subsys_anno, column_split = chrt$Clusters, row_split = subsys$SubSystem,
        heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                    title_gp = gpar(fontsize = 8), at = c(-1000, -10, 0, 10, 1000), labels_gp = gpar(fontsize = 8)), row_title = " ")

dev.off()


################################################################
