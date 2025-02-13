

# Load required packages
library(tidyverse)
library(rtables)
library(htmltools)

make_demographic_table <- function(metadata) {
  
  # Define race as a factor and arrange the values in the order needed in the table
  metadata <- metadata |> 
    mutate(race = factor(race)) |> 
    mutate(race = fct_relevel(race, c("white", "black or african american", "asian", "not reported")))
  metadata$race
  # Create the layout for demographic table
  # List of variables to be included in the rows
  vars <- c("age_at_index", "race")
  labels <- list("age_at_index" = "Age (years)",
                 "race" = "Race"
  )
  
  lyt <- basic_table() |>
    split_cols_by("gender") |>
    add_colcounts() |>
    analyze(vars, function(x) {
      if (is.numeric(x)) {
        in_rows(
          "Mean (SD)" = c(mean(x), sd(x)),
          "Median" = median(x),
          "Min - Max" = range(x),
          .formats = c("xx.x (xx.x)", "xx.x", "xx.x - xx.x")
        )
      } else if (is.factor(x) || is.character(x)) {
        in_rows(.list = list_wrap_x(table)(x))
      }
    },
    var_labels = labels)
  
  tbl <- build_table(lyt, metadata)
  #
  tbl_html <- as_html(tbl)
  #save_html(tbl_html, file = "PDAC-app/documents/demographic_table.html")
}
