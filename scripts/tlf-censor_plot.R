library(tidyverse)
library(cowplot)


make_censor_plot <- function(metadata) {
  # Define race as a factor and arrange the values in the order needed in the table
  metadata <- metadata |> 
    mutate(race = factor(race)) |> 
    mutate(race = fct_relevel(race, c("white", "black or african american", "asian", "not reported")))
  
  # Create a metadata subset for survival analysis
  metadata_survival_subset <- metadata |>
    select(c(sample_id, gender, race, ethnicity, ABSOLUTE.Purity, RPPA.Clusters..76.High.Purity.Samples.Only., Follow.up.vital.status, Follow.up.days, Days.to.death, Censored.1.yes.0.no))
  
  # There are some NaN values in the censor column that needs to be collected
  metadata_survival_subset <- metadata_survival_subset |>
    mutate(Censored = ifelse(Follow.up.vital.status == "Dead", 0, 1)) |>
    mutate(Follow.up.days = ifelse(Follow.up.vital.status == "Dead", Days.to.death, Follow.up.days))
  
  
  # Reference: https://ocbe-uio.github.io/survomics/survomics.html#introduction
  
  # frequency table of the patients w.r.t. status, gender and ethnicity
  metadata_survival_subset %>%
    count(Follow.up.vital.status, gender, race) %>%
    group_by(Follow.up.vital.status) %>%
    mutate(percent = prop.table(n)*100)
  
  # censoring plot by gender
  ID <- 1:nrow(metadata_survival_subset)
  
  metadata_survival_subset |> 
    arrange(gender, Follow.up.vital.status, desc(Follow.up.days)) |>
    ggplot(
      aes(y = ID, x = Follow.up.days, colour = gender, shape = factor(Follow.up.vital.status))
    ) +
    #geom_segment(aes(x = Follow.up.days, y = ID, xend = 0, yend = ID)) +
    geom_point() +
    ggtitle("") +
    labs(x = "Days", y = "Patients") +
    scale_shape_discrete(name = "Status", labels = c("Censored", "Dead")) +
    scale_color_discrete(
      name = "Gender",
      labels = c("female", "male")
    ) +
    theme_cowplot() +
    #theme(legend.position = "top", legend.direction = "vertical") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      #axis.title.y = element_blank(),   # Remove y-axis label
      axis.text.y = element_blank(),    # Remove y-axis text
      axis.ticks.y = element_blank()   # Remove y-axis ticks
      #axis.line.y = element_blank()     # Remove y-axis line
    )
  
  
}
# Prepare the survival data for analysis


#metadata <- read.csv("PDAC-app/data/pdac_demographic_clinical_molecular_metadata.csv")

#ggsave("PDAC-app/censor_plot.tiff")

