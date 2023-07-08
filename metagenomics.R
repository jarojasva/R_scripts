# 1 - Can you compare most changed strains in the S1 and S3 timepoints in responders in a heat map? You need to analyze the metagenomics data for this.

library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(vegan)
library(ecodist)
library(limma)
library(matrixStats)

# Read the abundance table
abundance_table <- read.table("merged_abundance_table_species_bacerias.txt", sep = "\t", header = TRUE, row.names = 1)
respodS1S3 <- abundance_table %>% select(starts_with("R_"))
respodS1S3.2 <- respodS1S3 %>% select(matches("_S1_unmapped$|_S3_unmapped$"))
zero_var_rows <- which(rowVars(as.matrix(respodS1S3.2)) == 0)
respodS1S3.2_filtered <- respodS1S3.2[-zero_var_rows, ]

#normalized_data <- log2(respodS1S3.2 + 1)  # Apply a log transformation to the abundance values
metadata <- read.table("sampledata_S1_S3.txt", sep = "\t", header = TRUE)

# Create a numeric design matrix
design_matrix <- model.matrix(~ 0 + metadata$Time)
colnames(design_matrix) <- c("S1", "S3")

# Perform differential abundance analysis
fit <- lmFit(respodS1S3.2_filtered, design_matrix)
contrast_matrix <- makeContrasts(S3 - S1, levels = colnames(design_matrix))
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

# Get the top differentially abundant strains
# top_strains <- topTable(fit_ebayes, coef = 1, adjust.method = "fdr", sort.by = "none", number = Inf)
# top_strain_names <- rownames(top_strains)
# top_strain_data <- respodS1S3.2_filtered[top_strain_names, ]
# top_strain_data <- data.matrix(top_strain_data)
# heatmap(top_strain_data, Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "blue"))(100), scale = "row")

# Get the top differentially abundant strains with P.Value < 0.05
top_strains <- topTable(fit_ebayes, coef = 1, adjust.method = "fdr", sort.by = "none", number = Inf)
filtered_strains <- top_strains[top_strains$P.Value < 0.05, ]
top_strain_names_filtered <- rownames(filtered_strains)
top_strain_data_filtered <- respodS1S3.2_filtered[top_strain_names_filtered, ]
top_strain_data_filtered <- data.matrix(top_strain_data_filtered)
Heatmap(top_strain_data_filtered, cluster_columns = FALSE, row_names_gp = grid::gpar(fontsize = 10), col = colorRampPalette(c("white", "blue"))(100), name = "Normalized Relative Abundance")



# 2 - Also can you compare the strain engraftment rate between R and NR for all time points?


library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(vegan)
library(ecodist)
library(limma)
library(matrixStats)
library(ggplot2)
library(reshape2)

# Read the abundance table
abundance_table <- read.table("merged_abundance_table_species_bacerias.txt", sep = "\t", header = TRUE, row.names = 1)

# Separate the R and NR samples and eliminate taxa rows with zero values
R_NR <- abundance_table %>% select(starts_with("R_"), starts_with("NR_"))
zero_var_rows <- which(rowVars(as.matrix(R_NR)) == 0)
R_NR_filtered <- R_NR[-zero_var_rows, ]

# Filter the columns starting with "R_" and "NR_"
R_columns <- grep("^R_", colnames(R_NR_filtered))
NR_columns <- grep("^NR_", colnames(R_NR_filtered))

# Calculate strain engraftment rate for R samples
engraftment_R <- apply(R_NR_filtered[, R_columns], 1, function(x) sum(x > 0) / length(x))

# Calculate strain engraftment rate for NR samples
engraftment_NR <- apply(R_NR_filtered[, NR_columns], 1, function(x) sum(x > 0) / length(x))

# Create a data frame with the strain names and engraftment rates
engraftment_data <- data.frame(Strain = rownames(R_NR_filtered), Engraftment_R = engraftment_R, Engraftment_NR = engraftment_NR)

# Perform t-test to compare engraftment rates between R and NR
engraftment_p_values <- apply(R_NR_filtered, 1, function(x) {
  t.test(x[R_columns], x[NR_columns])$p.value
})

# Add the p-values to the engraftment data frame
engraftment_data$p_value <- engraftment_p_values

# Sort the data frame by p-value in ascending order
engraftment_data <- engraftment_data[order(engraftment_data$p_value), ]

# Filter engraftment rates with p_values < 0.05
filtered_engraftment_data <- engraftment_data[engraftment_data$p_value < 0.05, ]

# Remove the "p_value" column from the dataframe
filtered_engraftment_data2 <- select(filtered_engraftment_data, -p_value)
fed <- melt(filtered_engraftment_data2)
colnames(fed) <- c("Strain", "Type", "Engraftment")

# Create a bar plot with customized colors and rotated x-axis labels
engraftment_plot <- ggplot(fed, aes(x = Strain, y = Engraftment, fill = Type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c("Engraftment_R" = "blue", "Engraftment_NR" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14)) +
  labs(x = "Strain", y = "Engraftment Rate", title = "Engraftment Rates between R and NR")

# Display the plot
engraftment_plot








# 3 - Do a bray-Curtis dissimilarity analysis to donor fold-change between R and NR for all time points

library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(vegan)
library(ecodist)
library(limma)
library(matrixStats)

# Read the abundance table
abundance_table <- read.table("merged_abundance_table_species_bacerias.txt", sep = "\t", header = TRUE, row.names = 1)

#ANALYSIS R vs NR , Bray Curtis logFC > 1.5

R_NR <- abundance_table %>% select(starts_with("R_"), starts_with("NR_"))
zero_var_rows <- which(rowVars(as.matrix(R_NR)) == 0)
R_NR_filtered <- R_NR[-zero_var_rows, ]

#normalized_data <- log2(respodS1S3.2 + 1)  # Apply a log transformation to the abundance values
metadata_R_NR <- read.table("sampledata_R_NR.txt", sep = "\t", header = TRUE)

# Create a numeric design matrix
design_matrix_R_NR <- model.matrix(~ 0 + metadata_R_NR$Time)
colnames(design_matrix_R_NR) <- c("S1", "S2", "S3","S4")

# Perform differential abundance analysis
fit_R_NR <- lmFit(R_NR_filtered, design_matrix_R_NR)
contrast_matrix_R_NR <- makeContrasts(S3 - S1, levels = colnames(design_matrix_R_NR))
fit_contrast_R_NR <- contrasts.fit(fit_R_NR, contrast_matrix_R_NR)
fit_ebayes_R_NR <- eBayes(fit_contrast_R_NR)

# Con top_strains P.Value < 0.05
top_strains_R_NR <- topTable(fit_ebayes_R_NR, coef = 1, adjust.method = "fdr", sort.by = "none", number = Inf)
filtered_strains_R_NR <- top_strains_R_NR[top_strains_R_NR$logFC > 1.5, ]
top_strain_names_filtered_R_NR <- rownames(filtered_strains_R_NR)
top_strain_data_filtered_R_NR <- R_NR_filtered[top_strain_names_filtered_R_NR, ]
top_strain_data_filtered_R_NR <- data.matrix(top_strain_data_filtered_R_NR)
# Heatmap(top_strain_data_filtered_R_NR, cluster_columns = FALSE, row_names_gp = grid::gpar(fontsize = 10))

# Calcular la disimilitud de Bray-Curtis para cada punto de tiempo
bray_curtis_dist <- vegan::vegdist(t(top_strain_data_filtered_R_NR), method = "bray")

bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
dim(bray_curtis_pcoa_df)
dim(metadata_R_NR)

bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df, metadata_R_NR)
bray_curtis_plot_R_NR <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df, aes(x = pcoa1, y = pcoa2, color = Time, shape = Responder), size = 5) +
  labs(x = "PC1",
       y = "PC2",
       title = "Bray-Curtis PCoA Responder_NonResponder") + 
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 12))
bray_curtis_plot_R_NR


#ANALYSIS D vs NR , Bray Curtis logFC > 1.5

D_NR <- abundance_table %>% select(starts_with("D"), starts_with("NR_"))
zero_var_rows <- which(rowVars(as.matrix(D_NR)) == 0)
D_NR_filtered <- D_NR[-zero_var_rows, ]

#normalized_data <- log2(respodS1S3.2 + 1)  # Apply a log transformation to the abundance values
metadata_D_NR <- read.table("sampledata_D_NR.txt", sep = "\t", header = TRUE)

# Create a numeric design matrix
design_matrix_D_NR <- model.matrix(~ 0 + metadata_D_NR$Time)
colnames(design_matrix_D_NR) <- c("HD", "S1", "S2", "S3","S4")

# Perform differential abundance analysis
fit_D_NR <- lmFit(D_NR_filtered, design_matrix_D_NR)
contrast_matrix_D_NR <- makeContrasts(S3 - S1, levels = colnames(design_matrix_D_NR))
fit_contrast_D_NR <- contrasts.fit(fit_D_NR, contrast_matrix_D_NR)
fit_ebayes_D_NR <- eBayes(fit_contrast_D_NR)

# Con top_strains P.Value < 0.05
top_strains_D_NR <- topTable(fit_ebayes_D_NR, coef = 1, adjust.method = "fdr", sort.by = "none", number = Inf)
filtered_strains_D_NR <- top_strains_D_NR[top_strains_D_NR$logFC > 1.5, ]
top_strain_names_filtered_D_NR <- rownames(filtered_strains_D_NR)
top_strain_data_filtered_D_NR <- D_NR_filtered[top_strain_names_filtered_D_NR, ]
top_strain_data_filtered_D_NR <- data.matrix(top_strain_data_filtered_D_NR)
# Heatmap(top_strain_data_filtered_D_NR, cluster_columns = FALSE, row_names_gp = grid::gpar(fontsize = 10))

# Calcular la disimilitud de Bray-Curtis para cada punto de tiempo
bray_curtis_dist <- vegan::vegdist(t(top_strain_data_filtered_D_NR), method = "bray")

bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
dim(bray_curtis_pcoa_df)
dim(metadata_D_NR)

bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df, metadata_D_NR)
bray_curtis_plot_D_NR <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df, aes(x = pcoa1, y = pcoa2, color = Time, shape = Type), size = 5) +
  labs(x = "PC1",
       y = "PC2",
       title = "Bray-Curtis PCoA Donor_NonResponder") + 
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 12))
bray_curtis_plot_D_NR



#ANALYSIS D vs R , Bray Curtis logFC > 1.5

D_R <- abundance_table %>% select(starts_with("D"), starts_with("R_"))
zero_var_rows <- which(rowVars(as.matrix(D_R)) == 0)
D_R_filtered <- D_R[-zero_var_rows, ]

#normalized_data <- log2(respodS1S3.2 + 1)  # Apply a log transformation to the abundance values
metadata_D_R <- read.table("sampledata_D_R.txt", sep = "\t", header = TRUE)

# Create a numeric design matrix
design_matrix_D_R <- model.matrix(~ 0 + metadata_D_R$Time)
colnames(design_matrix_D_R) <- c("HD", "S1", "S2", "S3","S4")

# Perform differential abundance analysis
fit_D_R <- lmFit(D_R_filtered, design_matrix_D_R)
contrast_matrix_D_R <- makeContrasts(S3 - S1, levels = colnames(design_matrix_D_R))
fit_contrast_D_R <- contrasts.fit(fit_D_R, contrast_matrix_D_R)
fit_ebayes_D_R <- eBayes(fit_contrast_D_R)

# Con top_strains P.Value < 0.05
top_strains_D_R <- topTable(fit_ebayes_D_R, coef = 1, adjust.method = "fdr", sort.by = "none", number = Inf)
filtered_strains_D_R <- top_strains_D_R[top_strains_D_R$logFC > 1.5, ]
top_strain_names_filtered_D_R <- rownames(filtered_strains_D_R)
top_strain_data_filtered_D_R <- D_R_filtered[top_strain_names_filtered_D_R, ]
top_strain_data_filtered_D_R <- data.matrix(top_strain_data_filtered_D_R)
# Heatmap(top_strain_data_filtered_D_R, cluster_columns = FALSE, row_names_gp = grid::gpar(fontsize = 10))

# Calcular la disimilitud de Bray-Curtis para cada punto de tiempo
bray_curtis_dist <- vegan::vegdist(t(top_strain_data_filtered_D_R), method = "bray")

bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
dim(bray_curtis_pcoa_df)
dim(metadata_D_R)

bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df, metadata_D_R)
bray_curtis_plot_D_R <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df, aes(x = pcoa1, y = pcoa2, color = Time, shape = Type), size = 5) +
  labs(x = "PC1",
       y = "PC2",
       title = "Bray-Curtis PCoA Donor_Responder") + 
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 12))
bray_curtis_plot_D_R