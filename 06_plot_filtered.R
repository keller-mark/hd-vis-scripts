library(ggplot2)
library(tidyr)


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

