### Filtering the clinical data by removing columns having more than 20% NA's

# read in your data
clinical_data = read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/all_data.txt")

# calculate the percentage of missing values in each column
missing_percent = colSums(is.na(clinical_data)) / nrow(clinical_data) * 100

# create a vector of column names to remove columns having more than 20% NA's
cols_to_remove = names(clinical_data)[missing_percent >= 20]

# remove the columns from your data
clinical_data = clinical_data[, !(names(clinical_data) %in% cols_to_remove)]

write.table((clinical_data),file="/home/shweta/Desktop/ProjectChildren/Clinical Data/Filtered_clinical_data.xlsx", quote = FALSE,col.names = NA, sep = "\t")



#####################################################

### Filtering the clinical data by removing columns having more than 20% NA's

# read in your data
clinical_data = read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Clinical_continuous_data.txt")

# calculate the percentage of missing values in each column
missing_percent = colSums(is.na(clinical_data)) / nrow(clinical_data) * 100

# create a vector of column names to remove columns having more than 20% NA's
cols_to_remove = names(clinical_data)[missing_percent >= 20]

# remove the columns from your data
clinical_data = clinical_data[, !(names(clinical_data) %in% cols_to_remove)]

write.table((clinical_data),file="/home/shweta/Desktop/ProjectChildren/Clinical Data/filtered_continuous_data.txt", quote = FALSE,col.names = NA, sep = "\t")


########################################################


### Filtering the clinical data by removing columns having more than 20% NA's

# read in your data
categorical_data = read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/categorical_data.txt")

# calculate the percentage of missing values in each column
missing_percent = colSums(is.na(categorical_data)) / nrow(categorical_data) * 100

# create a vector of column names to remove
cols_to_remove = names(categorical_data)[missing_percent >= 20]

# remove the columns from your data
categorical_data = categorical_data[, !(names(categorical_data) %in% cols_to_remove)]

write.table((categorical_data),file="/home/shweta/Desktop/ProjectChildren/Clinical Data/Filtered_categorical_data.txt", quote = FALSE,col.names = NA, sep = "\t")


########################################################

### To check if the Data is normalized

conti_data=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Clinical_continuous_data.txt",header=TRUE,row.names = 1)

# Create a text file for the output
output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/Normal_distribution_for_continuous_data/outliers.txt", "w")

# Loop for all columns containing continuous variables
for (i in 2:42) {
  # Create a PDF file with a name corresponding to the variable
  pdf(paste0("/home/shweta/Desktop/ProjectChildren/Clinical Data/Normal_distribution_for_continuous_data/variable_", colnames(conti_data)[i], ".pdf"))
  
  # Create a histogram
  hist(conti_data[, i], main = paste("Histogram for", colnames(conti_data)[i]), xlab = "Value", ylab = "Frequency")
  
  # Create a boxplot
  boxplot(conti_data[, i], main = paste("Boxplot for", colnames(conti_data)[i]), ylab = "Value")
  
  # Create a QQ plot
  qqnorm(conti_data[,i], main=paste("QQ plot for", colnames(conti_data)[i]))
  qqline(conti_data[,i])
  
  # Save and close the PDF file
  dev.off()
  
  # Calculate z-scores
  z_scores = scale(conti_data[, i])
  
  # Write z-scores to a text file
  write.table(z_scores, file = paste0("/home/shweta/Desktop/ProjectChildren/Clinical Data/Z_score/z_scores_", colnames(conti_data)[i], ".txt"))
  
  # Identify outliers
  outliers = which(abs(z_scores) > 3)
  
  # Print the number of outliers
  out = paste("There are", length(outliers), "outliers in", colnames(conti_data)[i], "\n")
  
  # Write outliers to a text file
  cat(out, file = output_file)
}

# Close the output file
close(output_file)


################################################



### Kruskal Wallis Test

df3=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/filtered_continuous_data.txt")

# Define the column names
col_names <- names(df3)[2:27]

# Open a file for writing results
result_file <- file("/home/shweta/Desktop/ProjectChildren/Clinical Data/kruskal_results2.txt", "w")

# Loop through each column name in a vector
for (col_name in col_names) {
  # Perform the Kruskal-Wallis test for the current column
  kw_result <- kruskal.test(as.formula(paste0(col_name, " ~ Clusters")), data = df3)
  
  # get the p-values
  p.values <- kw_result$p.value
  
  # apply BH correction
  p.adj <- p.adjust(p.values, method = "BH")
  
  # get the absolute p-values that passed the BH correction
  p.abs.adj <- abs(p.values[p.adj <= 0.05])
  
  
  # Write the test results to the file
  writeLines(capture.output(kw_result, p.values, p.adj, p.abs.adj), result_file)
  writeLines("\n", result_file)
}

# Close the file
close(result_file)


# Open a file for writing results
median_result_file <- file("/home/shweta/Desktop/ProjectChildren/Clinical Data/kw_median_IQR_result.txt", "w")

for (col_name in col_names) {
  
  # Calculate medians for each group
  #medians <- aggregate(as.formula(paste0(col_name, " ~ Clusters")), data = df3, median)
  
  medians_iqrs = aggregate(df3[col_names], by = list(Clusters = df3$Clusters),
                           FUN = function(x) c(median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE)))
  
  
  # Write the test results to the file
  writeLines(capture.output(medians_iqrs), median_result_file)
  writeLines("\n", median_result_file)
}

# Close the file
close(median_result_file)


########################################################

### Mann Whitney U-test

result_file <- file("/home/shweta/Desktop/ProjectChildren/Clinical Data/mann_w_result.txt", "w")


for (i in 2:ncol(df3)) {
  p_values <- pairwise.wilcox.test(df3[, i], df3$Clusters, p.adjust.method = "BH")
  
  col_n <- paste("Results for variable:", names(df3)[i])
  
  writeLines(col_n, result_file)
  writeLines(capture.output(p_values), result_file)
  writeLines("\n", result_file)
}
close(result_file)


###########################################


### Chi square test

categoical_data=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Filtered_categorical_data.txt",header=TRUE)
View(categoical_data)

# Check the number of unique values in the column
predictor_cols = colnames(categoical_data)[-1]

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_results2.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categoical_data$Clusters, categoical_data[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  
  # get the p-values
  p.values <- chisq_result$p.value
  
  # apply BH correction
  p.adj <- p.adjust(p.values, method = "BH")
  
  # get the absolute p-values that passed the BH correction
  p.abs.adj <- abs(p.values[p.adj <= 0.05])
  
  writeLines(capture.output(result_string, p.values, p.adj, p.abs.adj), output_file)
  writeLines("\n", output_file)
}

# Close the output file
close(output_file)




###########################################


# Chi square test Group1

categorical_data_Group1=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group1(Cluster1_2).txt",header=TRUE)
View(categorical_data_Group1)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group1)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group1.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group1$Clusters, categorical_data_Group1[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)


#####################################################


# Chi square test Group2

categorical_data_Group2=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group2(Cluster1_3).txt",header=TRUE)
View(categorical_data_Group2)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group2)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group2.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group2$Clusters, categorical_data_Group2[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)


#####################################################


# Chi square test Group3

categorical_data_Group3=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group3(Cluster1_4).txt",header=TRUE)
View(categorical_data_Group3)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group3)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group3.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group3$Clusters, categorical_data_Group3[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)


#####################################################

# Chi square test Group4

categorical_data_Group4=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group4(Cluster2_3).txt",header=TRUE)
View(categorical_data_Group4)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group4)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group4.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group4$Clusters, categorical_data_Group4[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)


#####################################################


# Chi square test Group5

categorical_data_Group5=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group5(Cluster2_4).txt",header=TRUE)
View(categorical_data_Group5)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group5)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group5.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group5$Clusters, categorical_data_Group5[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)

#####################################################

# Chi square test Group6

categorical_data_Group6=read.delim("/home/shweta/Desktop/ProjectChildren/Clinical Data/Group6(Cluster3_4).txt",header=TRUE)
View(categorical_data_Group6)

# Check the number of unique values in the column
predictor_cols = colnames(categorical_data_Group6)[-1]

# write.table(chisq_result,"P_values_chi_Square test.txt")

output_file = file("/home/shweta/Desktop/ProjectChildren/Clinical Data/chi_square_group6.txt", "w")

# Loop over the predictor variables and perform chi-square tests
for (col in predictor_cols) {
  contingency_table = table(categorical_data_Group6$Clusters, categorical_data_Group6[[col]])
  chisq_result = chisq.test(contingency_table)
  result_string = paste("Chi-square test for", col, "p-value:", chisq_result$p.value, "\n")
  writeLines(result_string, output_file)
}

# Close the output file
close(output_file)

#####################################################

set.seed(123)

# load the dplyr package
library(dplyr)

# gather all the categorical columns into key-value pairs
df = categoical_data %>% 
  tidyr::gather(key = Column, value = Value, -Clusters)

# group the data frame by clusters and column name, and count the number of each category
count_df = df %>%
  group_by(Clusters, Column, Value) %>%
  summarise(count = n()) %>%
  ungroup()

# calculate the percentage of each category for each cluster
percent_df = count_df %>%
  group_by(Clusters, Column) %>%
  mutate(percent = count/sum(count) * 100) %>%
  arrange(Column)

# view the result
View(percent_df)

write.table((percent_df),file = "/home/shweta/Desktop/ProjectChildren/Clinical Data/percentage_categorical_data.txt", quote = FALSE,col.names = NA)


#########################################################
