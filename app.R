# Import required packages
library(shiny)
library(bslib)
library(tidyverse)
library(rtables)
library(plotly)
library(DT)
library(shinyWidgets)

# Save today's date as a character string to use later in the name of the app
today <- format(Sys.Date(), format = "%d %B %Y")

# Load the data - we only need it once, so keeping it outside the ui and server functions
# The path of the files had to be given relative to working directory and not relative to the directory where app.R is located.
metadata <- read.csv("./data/pdac_demographic_clinical_molecular_metadata.csv")

rppa <- read.csv("./data/pdac_rppa.csv")

# Deal with NA values in protein data
# Remove rows with missing data for more than 10 samples
rppa <- rppa[rowSums(is.na(rppa))<10,]
rownames(rppa) <- rppa[,1]

proteomics <- data.frame(t(rppa[,6:103]))


# Calculate column means as a list
#Fill in missing values based on colMeans -Mean of that protein across all biological samples

col_means_list <- as.list(colMeans(proteomics, na.rm = TRUE))
proteomics <- replace_na(proteomics, col_means_list)

# Source the helper file for demographic table
source(file = "./scripts/tlf-demographic_table.R")

# source the helper file for censor plot
source(file = "./scripts/tlf-censor_plot.R")

# source the helper file for survival plot
source(file = "./scripts/tlf-survival_plot.R")
# Define the available variables to compare survival probabilities

# Convert all character columns to factors
#mutate(across(c(A, B, C, E), as.factor)) # select columns A to C, and E (by name)
metadata <- metadata |> mutate(across(where(is.character), as.factor))
metadata <- metadata |> mutate(KRAS.Mutated..1.or.0. = as.factor(KRAS.Mutated..1.or.0.))
survival_choices <- c("Tissue Purity" = "Purity.Class..high.or.low.", 
                      "KRAS mutation status" = "KRAS.Mutated..1.or.0.",
                      "Clinical Pathology grade" = "Clinical.pathologic.M", 
                      "History of chronic pancreatitis"="History.of.chronic.pancreatitis",
                      "Race" = "race", 
                      "Gender" = "gender",
                      "Protein Cluster" = "RPPA.Clusters..76.High.Purity.Samples.Only.")

# Source the helper file for PCA plot
source(file = "./scripts/tlf-pca-plot.R")

# Source the helper file for purity plots
source(file = "./scripts/tlf-clinical_overview.R")


# Setting up the UI
ui <- fluidPage(
  fluidRow(
    column(width = 12, align = "center", h1("Characterization of Pancreatic Ductal Adenocarcinoma (PDAC) data from TCGA"))
  ),
  fluidRow(
    column(width = 12, align = "center", p(today))
  ),
  page_navbar(
    title = NULL,
    nav_panel(
      title = "Information",
      includeMarkdown("./documents/information.md")
    ),
    nav_panel(
      title = "Demographic Table",
      h3("Summary of Demographic Characteristics", align = "center"),
      column(width = 6, offset = 3, uiOutput(outputId = "demographic_table"))
    ),
    nav_panel(
      title = "Censoring Plot",
      h3("Censoring plot of the subjects", align = "center"),
      column(width = 6, offset = 3,plotOutput(outputId = "censor_plot")),
      tags$footer(textOutput(outputId = "censor_plot_note"))
    ),
    nav_panel(
      title = "Kaplan Meier plot for survival",
      layout_sidebar(
        # insert a sidebar for filtering the data
        sidebar = sidebar(
          selectInput(inputId = "select_factor", 
                      label = "Select a variable to compare survival probabilities", 
                      choices = c("None",survival_choices), selected = "None")
        ),
      h3("Survival Plot and Risk Table", align = "center"),
      column(width = 6, offset = 3, plotOutput(outputId = "survival_plot")),
      column(width = 12, align = "center", textOutput(outputId = "log_rank_test_title")),
      # Add a table showing p-values of tests between different curves.
      column(width = 12, align = "center", tableOutput(outputId = "survival_comparison_table"))
      #textOutput(outputId = "survival_text")
    )),
    nav_panel(
      title = "Principal Components Analysis",
      layout_sidebar(
        # Insert sidebar for coloring PCA plot
        sidebar = sidebar(
          selectInput(inputId = "select_pca_color",
                      label = "Select color for PCA plot",
                      choices = survival_choices,
                      selected = "gender")
        ),
        column(width = 8, offset = 2, plotOutput(outputId = "pca_plot"))
      )),
    nav_panel(
      title = "Tumor Purity Analysis",
        navset_tab(
          nav_panel(
            title = "Overlap between molecular clusters, sample purity, and leukocyte fraction",
            alert( "Note: The plot takes a few seconds to render.",status = "info",
                  dismissible = TRUE),
            column(width = 12, offset = 0, plotOutput(outputId = "radial_plot")),
            column(width = 8, offset =2, uiOutput(outputId = "text_radial"))
          ),
          nav_panel(
            title = "Sample purity for different molecular classifications",
            layout_column_wrap(
              card(plotOutput(outputId = "purity_bailey_clusters"),
                   textOutput(outputId = "text_bailey")),
              card(plotOutput(outputId = "purity_moffitt_clusters"),
                   textOutput(outputId = "text_moffitt")),
              card(plotOutput(outputId = "purity_collisson_clusters"),
                   textOutput(outputId = "text_collisson")),
              card(plotOutput(outputId = "estimated_tumor_fraction"),
                   textOutput(outputId = "text_tumor_fraction")),
              width = 1/2
            )
            
          )
      )
    )
  )
)



# Define server logic
server <- function(input, output) {
  
  # Including a reactive HTML element for demographic table in the UI
  # In this case, I don't need to make the table as interactive HTML element because it is not changing
  output$demographic_table <- renderUI(
    {make_demographic_table(metadata)}
  )
  
  # Output the censor plot. The censor plot ggplot object is defined in a separate helper file
  output$censor_plot <- renderPlot({make_censor_plot(metadata)})
  
  # Output the note to print on the censor plot
  output$censor_plot_note <- renderText({
    "Note: Two samples for which follow up time was not recorded have been excluded from the plot"
  })
  
  # Create a reactive variable based on the user selection
  user_selection <- reactive({
    if (input$select_factor == "None"){
      "1"
    } else {
      input$select_factor
    }
  })
 #output$survival_text <- renderText({
#    paste("User selected variable:", user_selection())
#  })
  # Output the survival plot. The survival plot ggplot object is defined in a separate file
 output$survival_plot <- renderPlot({
    make_survival_plot(metadata, sep_variable = user_selection())
  })
 
 output$survival_comparison_table <- renderTable({
   surv_diff_table <- make_survival_table(metadata, sep_variable = user_selection())
   if (length(surv_diff_table) != 0){
    surv_diff_table$p.value
   }
   else {
     data.frame()
   }
   
 }, rownames = TRUE, colnames = TRUE, align = "c")
 
 output$log_rank_test_title <- renderText({
   surv_diff_table <- make_survival_table(metadata, sep_variable = user_selection())
   if (length(surv_diff_table) != 0){
     paste("Pairwise",surv_diff_table$method, "p.value table with multiple testing correction using", surv_diff_table$p.adjust.method, "method")
   }
 })
 
 # Create the PCA plot with the data colored by user-selected variable 
 output$pca_plot <- renderPlot({
   create_pca_plot(proteomics, metadata, input$select_pca_color)
 })
 
 # Create plots on purity analysis page
 output$purity_bailey_clusters <- renderPlot({
   create_bailey_purity_boxplot(metadata)
 })
 
 output$purity_moffitt_clusters <- renderPlot({
   create_moffitt_boxplot(metadata)
 })
 output$purity_collisson_clusters <- renderPlot({
   create_collisson_boxplot(metadata)
 })
 output$estimated_tumor_fraction <- renderPlot({
   create_tumor_fraction_plot(metadata)
 })
 
 output$radial_plot <- renderPlot({
   create_radial_plot(metadata)
 })
 
 #output$radial_plot_text <- renderText({
#   "The figure shows subjects classiffied according to gene signatures specified by Bailey et al. (outer circle), Collisson et al., Moffitt et al., Tumor Purity classification and estimated leukocyte percentage (innermost circle)"
# })
 
 output$text_bailey <- renderText({
   "Classification based on gene expression signatures published by Bailey et al. was correlated with sample purity. Samples classified as ADEX or Immunogenic were generally lower in purity."
 })
 output$text_moffitt <- renderText({
   "Classification based on gene expression signatures published by Moffitt et al. was not correlated with sample purity."
 })
 output$text_collisson <- renderText({
   "Classification based on gene expression signatures published by Collisson et al. was correlated with sample purity. Samples classified as Quasi-mesencymal and Exocrine-like were lower in purity."
 })
 output$text_tumor_fraction <- renderText({
   "Plot shows estimates of tumor purity in each sample using pathology review data and ABSOLUTE algorithm. The samples were classified as high purity based on ABSOLUTE value > 33%."
 })
  
 output$text_radial <- renderUI({
   HTML("<p>Radial plot showing sample overlap of classification based on gene expression signatures published by Bailey et al., Collisson et al., Moffitt et al., Sample purity and estimated Leukocyte percentage (outer circle to inner circle, respectively).</p> 
   <li>The pancreatic progenitor subtye of Bailey et al. overlaps with classical subtypes from Moffitt et al. and Collisson et al. </li>
   <li>Squamous subtype of Bailey et al. overlaps with Basal-like subtype of Mofitt et al.</li>
    <li> High purity samples overlap with Progenitor and Squamous subtypes from Bailey et al. and Classical subtype from Collisson et al.</li>
    <li> Immunogenic subtype from bailey et al. tends to have higher estimated leukocyte percentage and lower purity.</li>"
     )
 })
 
  }




# Run the app ----------
shinyApp(ui = ui, server = server)