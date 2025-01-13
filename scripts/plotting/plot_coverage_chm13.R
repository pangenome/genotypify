library(ggplot2)
library(dplyr)
library(tidyr)

# Read the depth data
depth_data <- read.table(
  "/home/guarracino/Desktop/Allentoft2024.ERR12119480.mm2.map.depth.windows.10kbp.bed",
  header=FALSE,
  col.names=c("chrom", "start", "end", "depth"),
  comment.char = '?'
)

# Create chromosome order (matching your fai file)
chrom_order <- c(
  "chm13#1#chr1", "chm13#1#chr2", "chm13#1#chr3", "chm13#1#chr4", "chm13#1#chr5",
  "chm13#1#chr6", "chm13#1#chr7", "chm13#1#chr8", "chm13#1#chr9", "chm13#1#chr10",
  "chm13#1#chr11", "chm13#1#chr12", "chm13#1#chr13", "chm13#1#chr14", "chm13#1#chr15",
  "chm13#1#chr16", "chm13#1#chr17", "chm13#1#chr18", "chm13#1#chr19", "chm13#1#chr20",
  "chm13#1#chr21", "chm13#1#chr22", "chm13#1#chrX", "chm13#1#chrY", "chm13#1#chrM"
)

# Convert chromosome to factor with correct order
depth_data$chrom <- factor(depth_data$chrom, levels=chrom_order)

# Calculate chromosome midpoints for x-axis labels
chrom_sizes <- c(
  248387328, 242696752, 201105948, 193574945, 182045439,
  172126628, 160567428, 146259331, 150617247, 134758134,
  135127769, 133324548, 113566686, 101161492, 99753195,
  96330374, 84276897, 80542538, 61707364, 66210255,
  45090682, 51324926, 154259566, 62460029, 16569
)
names(chrom_sizes) <- chrom_order

depth_data <- depth_data %>% filter(chrom != 'chm13#1#chrM')

# Create the plot
p <- ggplot(depth_data, aes(x=start, y=depth)) +
  # Add depth line
  geom_line(size=0.5, alpha=0.7) +
  
  # Facet by chromosome
  facet_wrap(~chrom, ncol=3, strip.position="top") + #, scales="free_x") +
  
  # Customize theme
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=16),
    axis.text.x = element_text(angle=0, size=16),
    axis.text.y = element_text(size=16),
    axis.title = element_text(size=20)
  ) +
  
  # Labels
  labs(
    x = "Position (Mbp)",
    y = "Average depth",
    title = "Read alignment depth across chromosomes"
  ) +
  
  # Format x-axis to show positions in Mb
  scale_x_continuous(labels = function(x) paste0(x/1e6), 
                     breaks = scales::pretty_breaks(n=5))

# If you want to add a horizontal line at mean depth
mean_depth <- mean(depth_data$depth, na.rm=TRUE)
p <- p + geom_hline(yintercept=mean_depth, linetype="dashed", color="red", alpha=0.5)

# If you want to highlight regions with unusual depth
# For example, regions >2x mean depth
p <- p + geom_rect(data=subset(depth_data, depth > 2*mean_depth),
                   aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
                   fill="red", alpha=0.2)

# If you want to log-scale the y-axis (useful for very variable depth)
#p <- p + scale_y_log10()

# Save the plot
ggsave("depth_plot.pdf", p, width=70, height=50, limitsize = FALSE)
