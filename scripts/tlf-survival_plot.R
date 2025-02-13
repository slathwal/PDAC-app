library(tidyverse)
library(survival)
library(ggsurvfit)
library(cowplot)
library(survminer)

# [TO DO] Add a table showing the p-values and confidence intervals of the survival probabilities by factors
# Use the appropriate formulae based on the factor - if column is numeric, we will need cox-PH test
# if column is factor we can use pairwise log-rank test with multiple p-value correction

make_survival_plot <- function(metadata, sep_variable) {
  
  # Define race as a factor and arrange the values in the order needed in the table
  metadata <- metadata |> 
    mutate(race = factor(race)) |> 
    mutate(race = fct_relevel(race, c("white", "black or african american", "asian", "not reported")))
  
  # Create a metadata subset for survival analysis
  metadata_survival_subset <- metadata |>
    select(c(sample_id, gender, race, ethnicity, ABSOLUTE.Purity, RPPA.Clusters..76.High.Purity.Samples.Only., Follow.up.vital.status, Follow.up.days, Days.to.death, Censored.1.yes.0.no,
             Purity.Class..high.or.low.,KRAS.Mutated..1.or.0.,Clinical.pathologic.M,History.of.chronic.pancreatitis))
 
  # There are some NaN values in the censor column that needs to be corrected, create a new column called censored
  metadata_survival_subset <- metadata_survival_subset |>
    mutate(Censored = ifelse(Follow.up.vital.status == "Dead", 1, 0)) |>
    mutate(Follow.up.days, Follow.up.days = ifelse(Follow.up.vital.status == "Dead", Days.to.death, Follow.up.days))
    
  # Note: the Surv() function in the {survival} package accepts by default TRUE/FALSE, where TRUE is event and FALSE is censored; 1/0 where 1 is event and 0 is censored; or 2/1 where 2 is event and 1 is censored. Please take care to ensure the event indicator is properly formatted.
  
  # The Surv() function creates a survival object to use as response in a model formula
  # the survfit() function creates a survival curve using Kaplan Meier method
  
  # Create a curve for the entire dataset
  #s1 <- survfit(Surv(Follow.up.days, Censored) ~ 1, data = metadata_survival_subset)
  #str(s1)
  
  # use ggsurvfit package to create the KM survival plot
  # Create the formula to give to surv() function
  surv_formula <- as.formula(paste("Surv(Follow.up.days, Censored) ~", sep_variable))
  #print(surv_formula)
  survfit2(surv_formula, data = metadata_survival_subset) |>
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall Survival Probability"
    ) +
    #add_confidence_interval() +
    add_risktable(size = 4.5) +
    theme_cowplot()
}
#selected_variable <- "race"
#make_survival_plot(metadata, sep_variable = selected_variable)

#selected_variable <- "1"
#make_survival_plot(metadata, sep_variable = selected_variable)

make_survival_table <- function(metadata, sep_variable){
  # Define race as a factor and arrange the values in the order needed in the table
  metadata <- metadata |> 
    mutate(race = factor(race)) |> 
    mutate(race = fct_relevel(race, c("white", "black or african american", "asian", "not reported")))
  
  # Create a metadata subset for survival analysis
  metadata_survival_subset <- metadata |>
    select(c(sample_id, gender, race, ethnicity, ABSOLUTE.Purity, RPPA.Clusters..76.High.Purity.Samples.Only., Follow.up.vital.status, Follow.up.days, Days.to.death, Censored.1.yes.0.no,
             Purity.Class..high.or.low.,KRAS.Mutated..1.or.0.,Clinical.pathologic.M,History.of.chronic.pancreatitis, RPPA.Clusters..All.150.Samples.))
  
  # There are some NaN values in the censor column that needs to be corrected, create a new column called censored
  metadata_survival_subset <- metadata_survival_subset |>
    mutate(Censored = ifelse(Follow.up.vital.status == "Dead", 1, 0)) |>
    mutate(Follow.up.days, Follow.up.days = ifelse(Follow.up.vital.status == "Dead", Days.to.death, Follow.up.days))
  
  # Create a table for comparing survival differences
  # Create the formula to give to surv() function
  surv_formula <- as.formula(paste("Surv(Follow.up.days, Censored) ~", sep_variable))
  if (is.character(metadata_survival_subset[[sep_variable]]) | is.factor(metadata_survival_subset[[sep_variable]])){
    survminer::pairwise_survdiff(surv_formula, data = metadata_survival_subset)
  }else {
    list()
  }
  #survival::survdiff(surv_formula, data = metadata_survival_subset)
}
#surv_diff_table <- make_survival_table(metadata, sep_variable = "KRAS.Mutated..1.or.0.")
#surv_diff_table
#length(surv_diff_table)
#str(surv_diff_table)
#surv_diff_table$method
#surv_diff_table$data.name
#surv_diff_table$p.value
#surv_diff_table$p.adjust.method
