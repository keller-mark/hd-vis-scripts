library(ggplot2)
library(tidyr)
library(dplyr)


# Plot papers by field
by_field_df <- read.csv("filtered_papers_by_field.csv", row.names = 1)
by_field_df <- by_field_df[order(by_field_df$paper_count),]
by_field_df$field <- factor(by_field_df$field, levels = by_field_df$field)

q <- ggplot(by_field_df, aes(y = field, x = paper_count)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of papers citing selected dimensionality reduction methods, 2018-2023",
    x = "Number of Papers",
    y = "Field"
  )
ggsave("filtered_papers_by_field.png", plot = q, width = 10, height = 7, dpi = 300)


# Plot papers by field and cited method
by_field_method_df <- read.csv("filtered_papers_by_field_and_method.csv", row.names = 1)
#by_field_method_df <- by_field_method_df[by_field_method_df$field != "Computer Science", ]
by_field_method_df <- complete(by_field_method_df, method_acronym, field)
by_field_method_df <- arrange(by_field_method_df, field, desc(method_acronym), paper_count)
by_field_method_df$method_acronym <- factor(by_field_method_df$method_acronym, levels = unique(by_field_method_df$method_acronym))

q <- ggplot(by_field_method_df, aes(y = method_acronym, x = field, fill = paper_count)) +
  geom_tile() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient(low="#dddddd", high="black", limits = c(0, 1000), na.value = "white", oob = scales::squish) +
  labs(
    title = "Number of papers citing selected dimensionality reduction methods, 2018-2023",
    x = "Field",
    y = "Method"
  )

ggsave("filtered_papers_by_field_and_method.png", plot = q, width = 10, height = 15, dpi = 300)


# Plot papers by field and cited method, normalized within each field
norm_by_field_method_df <- (by_field_method_df
    %>% group_by(field)
    %>% mutate(norm_paper_count = paper_count/sum(paper_count, na.rm = TRUE))
   )

# Heatmap version
q <- ggplot(norm_by_field_method_df, aes(x = method_acronym, y = field, fill = norm_paper_count)) +
  geom_tile() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient(low="#dddddd", high="black", limits = c(0, 1), na.value = "white", oob = scales::squish) +
  labs(
    title = "Number of papers citing selected dimensionality reduction methods, 2018-2023",
    x = "Field",
    y = "Method"
  )

ggsave("filtered_papers_by_field_and_method_norm.png", plot = q, width = 15, height = 10, dpi = 300)

# Grouped bar plot (row per field) version
q <- ggplot(norm_by_field_method_df, aes(x = method_acronym, y = norm_paper_count)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(
    title = "Relative amount of papers by field citing selected dimensionality reduction methods, 2018-2023",
    x = "Method",
    y = "Normalized Paper Count"
  ) + facet_grid(rows = vars(field)) + theme(strip.text.y = element_text(angle = 0))

ggsave("filtered_papers_by_field_and_method_norm_bars.png", plot = q, width = 15, height = 10, dpi = 300)

# Stacked bar plot (stacked bar per field) version
max_method_norm_df <- norm_by_field_method_df %>% group_by(method_acronym) %>% summarise(max_norm_paper_count = max(norm_paper_count, na.rm = TRUE))
big_methods <- max_method_norm_df[max_method_norm_df$max_norm_paper_count >= 0.1, ]$method_acronym

big_norm_by_field_method_df <- norm_by_field_method_df[norm_by_field_method_df$method_acronym %in% big_methods, ]
q <- ggplot(big_norm_by_field_method_df, aes(x = field, y = norm_paper_count, fill = method_acronym)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(
    title = "Relative amount of papers by field citing selected dimensionality reduction methods, 2018-2023. Methods with < 10% usage omitted.",
    x = "Field",
    y = "Normalized Paper Count"
  )

ggsave("filtered_papers_by_field_and_method_norm_stacked_bar.png", plot = q, width = 15, height = 10, dpi = 300)



# Plot papers by venue
by_venue_df <- read.csv("filtered_papers_by_venue.csv", row.names = 1)
by_venue_df <- by_venue_df[order(by_venue_df$paper_count, decreasing = TRUE),]
by_venue_df$venue <- factor(by_venue_df$venue, levels = rev(by_venue_df$venue))
# Keep only top 50
by_venue_df <- by_venue_df[1:50, ]

q <- ggplot(by_venue_df, aes(y = venue, x = paper_count)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of papers citing selected dimensionality reduction methods, 2018-2023, top 50 venues",
    x = "Number of Papers",
    y = "Venue"
  )
ggsave("filtered_papers_by_venue.png", plot = q, width = 13, height = 30, dpi = 300)

# Plot papers by venue and cited method
by_field_method_df <- read.csv("filtered_papers_by_venue_and_method.csv", row.names = 1)
by_field_method_df <- complete(by_field_method_df, method_acronym, venue)
by_field_method_df <- arrange(by_field_method_df, venue, desc(method_acronym), paper_count)
by_field_method_df$method_acronym <- factor(by_field_method_df$method_acronym, levels = unique(by_field_method_df$method_acronym))
# Keep only top 50
by_field_method_df <- by_field_method_df[by_field_method_df$venue %in% by_venue_df$venue, ]

q <- ggplot(by_field_method_df, aes(y = method_acronym, x = venue, fill = paper_count)) +
  geom_tile() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient(low="#dddddd", high="black", limits = c(0, 100), na.value = "white", oob = scales::squish) +
  labs(
    title = "Number of papers citing selected dimensionality reduction methods, 2018-2023, top 50 venues",
    x = "Venue",
    y = "Method"
  )

ggsave("filtered_papers_by_venue_and_method.png", plot = q, width = 10, height = 15, dpi = 300)

