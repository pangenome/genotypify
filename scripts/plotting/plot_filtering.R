library(ggplot2)
library(dplyr)
library(tidyr)
library(scales) # For pretty_breaks
library(stringr) # For string manipulation

# Read the data
data <- read.delim("/home/guarracino/Desktop/filtering.tsv", header = TRUE, sep = "\t")

# Make sure we're working with a factor to maintain order
data$step <- factor(data$step, levels = c("RAW", "SPAN", "FLAGGER"))

# Create a wide format to calculate the differences
data_wide <- data %>%
  pivot_wider(names_from = step, values_from = num_samples) %>%
  mutate(
    RAW_minus_SPAN = RAW - SPAN,
    RAW_minus_FLAGGER = RAW - FLAGGER
  )

# Filter regions where the difference between RAW and SPAN is > 0
filtered_regions <- data_wide %>%
  filter(RAW_minus_SPAN > 100) %>%
  pull(region)

# Filter the original data to include only those regions
filtered_data <- data %>%
  filter(region %in% filtered_regions)

# Convert to wide format for sorting
filtered_wide <- filtered_data %>%
  pivot_wider(names_from = step, values_from = num_samples) %>%
  mutate(
    RAW_minus_SPAN = RAW - SPAN,
    RAW_minus_FLAGGER = RAW - FLAGGER
  )

# Sort by RAW - SPAN first, then by RAW - FLAGGER
sorted_data <- filtered_wide %>%
  arrange(desc(RAW_minus_SPAN), desc(RAW_minus_FLAGGER))

# Get the regions in the sorted order
sorted_regions <- sorted_data$region

# Filter and reorder the original data based on the sorted regions
filtered_sorted_data <- filtered_data %>%
  mutate(region = factor(region, levels = sorted_regions))

# Filter out RAW from your filtered_sorted_data
filtered_sorted_data_no_raw <- filtered_sorted_data %>%
  filter(step != "RAW")

# Extract everything from the 4th element onwards after splitting by '_' to use as x-tick labels
region_labels <- sapply(levels(filtered_sorted_data_no_raw$region), function(x) {
  parts <- str_split(x, "_")[[1]]
  if(length(parts) >= 4) {
    # Concatenate all parts from the 4th element to the end
    return(paste(parts[4:length(parts)], collapse = "_"))
  } else {
    return(x) # If there are fewer than 4 parts, use the original name
  }
})

# Create the modified plot with strictly positive y-axis
p1 <- ggplot(filtered_sorted_data_no_raw, aes(x = region, y = num_samples, fill = step)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 400, linetype = "dashed", color = "black") +
  # Add the label for RAW line
  annotate("text", x = 1, y = 400, label = "RAW = 400", vjust = -0.5, hjust = 0, color = "black") +
  theme_minimal() +
  labs(
    x = "Genomic Region",
    y = "Number of Samples",
    fill = "Processing Step"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
    #axis.text.x = element_blank(),  # This hides the x-axis labels
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"),
    legend.position = "top"
  ) +
  scale_fill_brewer(palette = "Set2") +
  # Use expand to prevent negative values on y-axis
  scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(n = 20),
    expand = expansion(mult = c(0, 0.05))  # This prevents negative values
  ) +
  # Use the extracted 4th elements as x-axis labels
  scale_x_discrete(labels = region_labels)

# Display plot
print(p1)

# Print the sorted data with differences for examination
print(sorted_data %>% select(region, RAW, SPAN, FLAGGER, RAW_minus_SPAN, RAW_minus_FLAGGER))

