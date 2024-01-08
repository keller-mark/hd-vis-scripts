library(UpSetR)


l1_df <- read.csv(file.path("data", "queries", "q06_get_l1_papers_fields_with_cp.csv"), row.names = 1)
l1_df <- l1_df[l1_df$source == "s2-fos-model", ]
l1_df <- l1_df[!(l1_df$field %in% c("Computer Science", "Engineering", "Mathematics")), ]


all_fields <- unique(l1_df$field)

list_input <- list()

for(field in all_fields) {
  list_input[[field]] <- l1_df[l1_df$field == field, ]$corpus_id
}

upset(fromList(list_input), order.by = "freq", empty.intersections = "on", nsets = 15)
