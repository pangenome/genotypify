# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)

# Read the TSV file
# Note: Replace 'your_file.tsv' with actual filename
data <- read_tsv("/home/guarracino/Desktop/ancient2406.depth.windows.200kbp.bed.gz",
                 col_names = c("chromosome", "start", "end", "coverage", 
                               "label", "dataset", "file"),
                 col_types = cols(
                   chromosome = col_character(),
                   start = col_double(),
                   end = col_double(),
                   coverage = col_double(),
                   label = col_factor(),
                   dataset = col_character(),
                   file = col_character()
                 ))

# Create chromosome order (matching your fai file)
chrom_order <- c(
  "chm13#1#chr1", "chm13#1#chr2", "chm13#1#chr3", "chm13#1#chr4", "chm13#1#chr5",
  "chm13#1#chr6", "chm13#1#chr7", "chm13#1#chr8", "chm13#1#chr9", "chm13#1#chr10",
  "chm13#1#chr11", "chm13#1#chr12", "chm13#1#chr13", "chm13#1#chr14", "chm13#1#chr15",
  "chm13#1#chr16", "chm13#1#chr17", "chm13#1#chr18", "chm13#1#chr19", "chm13#1#chr20",
  "chm13#1#chr21", "chm13#1#chr22", "chm13#1#chrX", "chm13#1#chrY", "chm13#1#chrM"
)

# Convert chromosome to factor with correct order
data$chromosome <- factor(data$chromosome, levels=chrom_order)



# Calculate number of unique files per dataset
dataset_counts <- data %>%
  group_by(dataset) %>%
  summarise(n_files = n_distinct(file)) %>%
  mutate(dataset_label = paste0(dataset, "\n(n=", n_files, ")"))

# Add this information to the main dataset
data <- data %>%
  left_join(dataset_counts, by = "dataset")



p <- ggplot(data, aes(x = label, y = coverage + 1, fill = label)) +
  geom_violin(alpha = 0.7, scale = "width") +  # Changed to violin plot
  # Add individual points with jitter for better visualization
  #geom_jitter(width = 0.2, alpha = 0.3) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)(c(1,1000))) +
  theme_minimal() +
  labs(
    title = "Coverage Distribution by Region Type",
    x = "Region",
    y = "Coverage + 1",
    fill = "Region Type"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
#p
ggsave(p, file="coverage_boxplot.png", width = 8, height = 6, dpi = 300, bg = "white")


# Create the plot with the new labels
p2 <- ggplot(data, aes(x = label, y = coverage + 1, fill = label)) +
  geom_violin(alpha = 0.7, scale = "width") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)(c(1,1000))) +
  facet_wrap(~dataset_label, scales = "free_y") +  # Use the new label with counts
  theme_minimal() +
  labs(
    title = "Coverage Distribution by Region Type and Dataset",
    x = "Region",
    y = "Coverage + 1",
    fill = "Region Type"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )
#p2
ggsave(p2, file="coverage_boxplot.by-dataset.png", width = 15, height = 9, dpi = 300, bg = "white")


p3 <- ggplot(data, aes(x = label, y = coverage + 1, fill = label)) +
  geom_violin(alpha = 0.7, scale = "width") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)(c(1,1000))) +
  facet_wrap(~chromosome, scales = "free", ncol = 7) +  # Changed faceting to chromosome
  theme_minimal() +
  labs(
    title = "Coverage Distribution by Region Type and Chromosome",  # Updated title
    x = "Region",
    y = "Coverage + 1",
    #fill = "Region Type"
  ) +
  #scale_fill_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )
#p3

# Save the plot
ggsave(p3, file="coverage_boxplot.by-chromosome.png", width = 19, height = 12, dpi = 300, bg = "white")

# If you want to add statistical summary
summary_stats <- data %>%
  group_by(label) %>%
  summarise(
    mean_coverage = mean(coverage),
    median_coverage = median(coverage),
    sd_coverage = sd(coverage),
    n = n()
  )

# Print summary statistics
print(summary_stats)

# Optionally save the plot
# ggsave("coverage_boxplot.pdf", width = 8, height = 6)