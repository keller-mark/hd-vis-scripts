library(UpSetR)


with_field_df <- read.csv("filtered_papers_with_field.csv", row.names = 1)

all_fields <- unique(with_field_df$field)

list_input <- list()

for(field in all_fields) {
  list_input[[field]] <- with_field_df[with_field_df$field == field, ]$paper_id
}

upset(fromList(list_input), order.by = "freq", empty.intersections = "on", nsets = 15)