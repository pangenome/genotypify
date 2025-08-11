# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(ggcorrplot)
library(viridis)
library(gridExtra)
library(scales)

# Load the data with fill=TRUE to handle the trailing tab
data <- read.table("/home/guarracino/Desktop/summary.tsv", 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE,
                   fill = TRUE)
data <- data[!grepl("^HARsv2", data$name), ]

# Remove the last column if it's empty (created by trailing tab)
if(all(is.na(data[,ncol(data)]))) {
  data <- data[,-ncol(data)]
}

# Check the structure
str(data)
dim(data)

# Remove rows with missing values in key numeric columns
data_clean <- data %>%
  filter(!is.na(num.haps) & num.haps > 0)

# CREATE PROPER CHROMOSOME ORDERING
# Define the correct chromosome order
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# Convert chrom column to factor with correct ordering
data_clean <- data_clean %>%
  mutate(chrom = factor(chrom, levels = chrom_order))

# Set theme for all plots with white background
theme_set(theme_minimal(base_size = 12) + 
            theme(plot.background = element_rect(fill = "white", color = NA),
                  panel.background = element_rect(fill = "white", color = NA)))

# Define QV colors with thresholds
qv_colors <- c(
  "high" = "#4CAF50",
  "mid" = "#FFC107",
  "low" = "#FF8C00",
  "verylow" = "#F44336"
)

qv_labels <- c(
  "high" = "High (>33)",
  "mid" = "Mid (>23, ≤33)",
  "low" = "Low (>17, ≤23)",
  "verylow" = "Very Low (≤17)"
)

# 1. SINGLE COLUMN PLOTS

# Plot 1: Distribution of number of haplotypes
p1 <- ggplot(data_clean, aes(x = num.haps)) +
  geom_histogram(binwidth = 10, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(title = "Distribution of Number of Haplotypes",
       x = "Number of Haplotypes",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Plot 2: Average haplotype length by chromosome (NOW PROPERLY SORTED)
p2 <- ggplot(data_clean, aes(x = chrom, y = avg.hap.len)) +
  geom_boxplot(fill = "lightcoral", alpha = 0.7) +
  labs(title = "Average Haplotype Length by Chromosome",
       x = "Chromosome",
       y = "Average Haplotype Length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Plot 3: Number of clusters distribution
p3 <- ggplot(data_clean, aes(x = num.clusters)) +
  geom_bar(fill = "darkgreen", alpha = 0.7) +
  labs(title = "Distribution of Number of Clusters",
       x = "Number of Clusters",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Plot 4: Quality value distribution (stacked bar) - UPDATED WITH NEW COLORS
qv_data <- data_clean %>%
  select(name, starts_with("perc.haps.QV")) %>%
  pivot_longer(cols = -name, names_to = "QV_type", values_to = "percentage") %>%
  mutate(QV_type = gsub("perc.haps.QV.", "", QV_type))

qv_summary <- qv_data %>% 
  group_by(QV_type) %>% 
  summarise(mean_perc = mean(percentage, na.rm = TRUE)) %>%
  mutate(QV_label = qv_labels[QV_type])

p4 <- ggplot(qv_summary, aes(x = reorder(QV_label, -mean_perc), y = mean_perc, fill = QV_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = qv_colors) +
  labs(title = "Average Quality Value Distribution",
       x = "Quality Value Category",
       y = "Average Percentage") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Save individual plots with white background
ggsave("num_haplotypes_dist.png", p1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("avg_hap_length_by_chrom.png", p2, width = 12, height = 6, dpi = 300, bg = "white")
ggsave("num_clusters_dist.png", p3, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("quality_value_dist.png", p4, width = 10, height = 6, dpi = 300, bg = "white")

# 2. MULTIPLE COLUMNS WITH CURVES

# Plot 5: Relationship between average haplotype length and Jaccard distance
p5 <- ggplot(data_clean, aes(x = avg.hap.len, y = avg.jaccard.dist)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(title = "Average Haplotype Length vs Average Jaccard Distance",
       x = "Average Haplotype Length",
       y = "Average Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Plot 6: Multiple length metrics comparison
length_data <- data_clean %>%
  select(name, min.hap.len, max.hap.len, avg.hap.len) %>%
  pivot_longer(cols = -name, names_to = "length_type", values_to = "length") %>%
  mutate(length_type = gsub(".hap.len", "", length_type))

p6 <- ggplot(length_data, aes(x = length_type, y = length, fill = length_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Comparison of Haplotype Length Metrics",
       x = "Length Type",
       y = "Length") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  scale_y_log10(labels = comma)

# Plot 7: Trends across genomic position (NOW WITH PROPERLY SORTED CHROMOSOMES)
p7 <- ggplot(data_clean, 
             aes(x = start/1000000, y = avg.jaccard.dist, color = chrom)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~chrom, scales = "free_x") +
  labs(title = "Jaccard Distance Trends Across Genomic Positions",
       x = "Genomic Position (Mb)",
       y = "Average Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  scale_x_continuous(labels = comma)

# Plot 8: QV categories vs number of clusters
qv_cluster_data <- data_clean %>%
  mutate(qv_high_category = cut(perc.haps.QV.high, 
                                breaks = c(-Inf, 50, 90, 95, 100),
                                labels = c("<50%", "50-90%", "90-95%", "95-100%")))

p8 <- ggplot(qv_cluster_data, aes(x = num.clusters, fill = qv_high_category)) +
  geom_histogram(binwidth = 2, position = "stack") +
  scale_fill_viridis_d() +
  labs(title = "Number of Clusters by QV High Percentage Categories",
       x = "Number of Clusters",
       y = "Count",
       fill = "QV High %") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Save multi-column plots with white background
ggsave("hap_length_vs_jaccard.png", p5, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("length_metrics_comparison.png", p6, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("genomic_position_trends.png", p7, width = 14, height = 8, dpi = 300, bg = "white")
ggsave("clusters_by_qv_categories.png", p8, width = 10, height = 6, dpi = 300, bg = "white")

# 3. CORRELATION HEATMAP - FIXED VERSION

# Select numeric columns for correlation
numeric_cols <- data_clean %>%
  select(num.haps,
         ref.length, avg.hap.len, 
         avg.jaccard.dist,
         num.clusters,
         perc.haps.QV.high, perc.haps.QV.verylow)

# Remove columns with zero variance or all NA values
numeric_cols_clean <- numeric_cols[, sapply(numeric_cols, function(x) {
  !all(is.na(x)) && var(x, na.rm = TRUE) != 0
})]

# Calculate correlation matrix
cor_matrix <- cor(numeric_cols_clean, use = "pairwise.complete.obs", method="spearman")

# Check for any remaining NA/Inf values and replace with 0
cor_matrix[is.na(cor_matrix)] <- 0
cor_matrix[is.infinite(cor_matrix)] <- 0

# Plot 9: Correlation heatmap using ggcorrplot
p9 <- ggcorrplot(cor_matrix, 
                 type = "lower",
                 lab = TRUE,
                 lab_size = 8,
                 colors = c("#6D9EC1", "white", "#E46726"),
                 title = "Spearman Correlation Heatmap",
                 ggtheme = theme_minimal()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 16),          # Increase legend text size
        legend.title = element_text(size = 20, face = "bold"),  # Increase legend title size
        legend.key.size = unit(1.5, "cm"),              # Increase legend key size
        legend.key.width = unit(1, "cm"),               # Adjust legend key width
        legend.key.height = unit(2, "cm"))              # Adjust legend key height
ggsave("correlation_heatmap.png", p9, width = 12, height = 10, dpi = 300, bg = "white")

# Alternative correlation heatmap with hierarchical clustering - FIXED
if(nrow(cor_matrix) > 1 && !any(is.na(cor_matrix))) {
  library(pheatmap)
  png("correlation_heatmap_clustered.png", width = 12, height = 10, units = "in", res = 300, bg = "white")
  pheatmap(cor_matrix,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           display_numbers = TRUE,
           number_format = "%.2f",
           number_color = "black",
           fontsize_number = 8,
           main = "Hierarchical Clustered Correlation Heatmap",
           clustering_method = "complete")
  dev.off()
}

# 4. ADDITIONAL ANALYSIS PLOTS

# Plot 10: Scatter plot matrix for key variables
library(GGally)
key_vars <- data_clean %>%
  select(avg.hap.len, sd.hap.length, avg.jaccard.dist, num.clusters, perc.haps.QV.high) %>%
  filter(complete.cases(.))  # Remove rows with NA values

if(nrow(key_vars) > 10) {  # Only create plot if we have enough data
  p10 <- ggpairs(key_vars,
                 upper = list(continuous = wrap("cor", size = 3)),
                 lower = list(continuous = wrap("points", alpha = 0.3)),
                 diag = list(continuous = wrap("densityDiag", fill = "lightblue")),
                 title = "Scatter Plot Matrix of Key Variables") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave("scatter_matrix.png", p10, width = 12, height = 10, dpi = 300, bg = "white")
}

# Plot 11: Time series-like plot for genomic regions (NOW PROPERLY SORTED)
p11 <- data_clean %>%
  arrange(chrom, start) %>%
  mutate(index = row_number()) %>%
  ggplot(aes(x = index, y = avg.jaccard.dist)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(color = chrom), size = 0.5) +
  labs(title = "Average Jaccard Distance Across All Genomic Regions",
       x = "Genomic Region Index (sorted by chromosome and position)",
       y = "Average Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  guides(color = guide_legend(nrow = 2))

ggsave("jaccard_distance_series.png", p11, width = 14, height = 6, dpi = 300, bg = "white")

# Plot 12: Faceted density plots for QV distributions - UPDATED WITH NEW COLORS
qv_density_data <- data_clean %>%
  select(starts_with("perc.haps.QV")) %>%
  pivot_longer(everything(), names_to = "QV_type", values_to = "percentage") %>%
  mutate(QV_type = gsub("perc.haps.QV.", "", QV_type),
         QV_label = qv_labels[QV_type])

p12 <- ggplot(qv_density_data, aes(x = percentage, fill = QV_type)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~QV_label, scales = "free_y") +
  scale_fill_manual(values = qv_colors) +
  labs(title = "Density Distributions of Quality Value Percentages",
       x = "Percentage",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("qv_density_distributions.png", p12, width = 12, height = 8, dpi = 300, bg = "white")

# Create a stacked bar chart showing QV distribution across all regions
p13 <- data_clean %>%
  select(name, starts_with("perc.haps.QV")) %>%
  pivot_longer(cols = -name, names_to = "QV_type", values_to = "percentage") %>%
  mutate(QV_type = gsub("perc.haps.QV.", "", QV_type),
         QV_label = factor(qv_labels[QV_type], levels = c("High (>33)", "Mid (>23, ≤33)", "Low (>17, ≤23)", "Very Low (≤17)"))) %>%
  group_by(name) %>%
  arrange(desc(QV_label)) %>%
  ggplot(aes(x = reorder(name, percentage), y = percentage, fill = QV_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = qv_colors, labels = qv_labels) +
  labs(title = "Quality Value Distribution Across All Genomic Regions",
       x = "Genomic Region",
       y = "Percentage",
       fill = "Quality Value") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  coord_cartesian(ylim = c(0, 100))

ggsave("qv_stacked_all_regions.png", p13, width = 14, height = 8, dpi = 300, bg = "white")

# Create a plot showing regions with problematic quality values
p14 <- data_clean %>%
  mutate(problematic = perc.haps.QV.low + perc.haps.QV.verylow > 20) %>%
  arrange(desc(perc.haps.QV.low + perc.haps.QV.verylow)) %>%
  head(50) %>%
  ggplot(aes(x = reorder(name, perc.haps.QV.low + perc.haps.QV.verylow), 
             y = perc.haps.QV.low + perc.haps.QV.verylow)) +
  geom_bar(stat = "identity", fill = "#F44336") +
  labs(title = "Top 50 Regions with Highest Low/Very Low Quality Values",
       x = "Genomic Region",
       y = "Combined Low + Very Low QV Percentage") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "black", alpha = 0.5)

ggsave("problematic_qv_regions.png", p14, width = 14, height = 8, dpi = 300, bg = "white")

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
print(summary(data_clean[, c("num.haps", "avg.hap.len", "avg.jaccard.dist", 
                             "num.clusters", "perc.haps.QV.high")]))

# Print correlation of key variables
cat("\n=== KEY CORRELATIONS ===\n")
key_vars_for_cor <- data_clean[, c("avg.hap.len", "avg.jaccard.dist", "num.clusters", "perc.haps.QV.high")]
if(all(sapply(key_vars_for_cor, function(x) var(x, na.rm = TRUE) > 0))) {
  cor_key <- cor(key_vars_for_cor, use = "pairwise.complete.obs")
  print(round(cor_key, 3))
}

# Print QV distribution summary
cat("\n=== QUALITY VALUE DISTRIBUTION SUMMARY ===\n")
qv_summary_stats <- data_clean %>%
  summarise(
    `High (>33) %` = mean(perc.haps.QV.high, na.rm = TRUE),
    `Mid (>23, ≤33) %` = mean(perc.haps.QV.mid, na.rm = TRUE),
    `Low (>17, ≤23) %` = mean(perc.haps.QV.low, na.rm = TRUE),
    `Very Low (≤17) %` = mean(perc.haps.QV.verylow, na.rm = TRUE)
  )
print(round(qv_summary_stats, 2))

# Identify problematic regions
cat("\n=== REGIONS WITH >20% LOW/VERY LOW QUALITY VALUES ===\n")
problematic_regions <- data_clean %>%
  mutate(low_quality_perc = perc.haps.QV.low + perc.haps.QV.verylow) %>%
  filter(low_quality_perc > 20) %>%
  select(name, chrom, start, end, low_quality_perc, num.haps, avg.jaccard.dist) %>%
  arrange(desc(low_quality_perc))
print(head(problematic_regions, 20))