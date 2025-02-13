library(ComplexHeatmap)
library(gridBase)
library(circlize)
library(ggplot2)
library(tidyr)
library(cowplot)


fix_data <- function(metadata){
  md <- metadata |> select(c(sample_id,Tumor.Sample.ID,ABSOLUTE.Purity,Purity.Class..high.or.low.,
                             Pathologist.Reviewed.Tumor.Cellularity,Initial.Slide.Tumor.Cellularity,
                             DNA.methylation.leukocyte.percent.estimate,mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical,
                             mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM,
                             mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX)) |>
    mutate(bailey_clusters = case_when(
      mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == 1 ~ "Squamous",
      mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == 2 ~ "Immunogenic",
      mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == 3 ~ "Progenitor",
      mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == 4 ~ "ADEX"
    ) ) |>
    mutate(moffit_clusters = case_when(
      mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical == 1 ~ "Basal-like",
      mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical == 2 ~ "Classical"
    )) |>
    mutate(collisson_clusters = case_when(
      mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM == 1 ~ "Classical",
      mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM == 2 ~ "Exocrine-like",
      mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM == 3 ~ "Quasi-mesenchymal"
    ))
  
}

create_bailey_purity_boxplot <- function(metadata){
  md <- fix_data(metadata)
  ggplot(data = md) +
    geom_boxplot(mapping= aes(x = bailey_clusters,
                              y = ABSOLUTE.Purity, colour = bailey_clusters),staplewidth = 0.5,
                 outliers = FALSE) +
    geom_point(mapping = aes(x = bailey_clusters, y = ABSOLUTE.Purity),
               position = position_jitter(width = 0.1, height = 0.1)) +
    theme_cowplot() +
    labs(x = NULL) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
}

create_moffitt_boxplot <- function(metadata){
  md <- fix_data(metadata)
  ggplot(data = md) +
    geom_boxplot(mapping= aes(x = moffit_clusters,
                              y = ABSOLUTE.Purity, colour = moffit_clusters), staplewidth = 0.5,
                 outliers = FALSE) +
    geom_point(mapping = aes(x = moffit_clusters, y = ABSOLUTE.Purity),
               position = position_jitter(width = 0.1, height = 0.1)) +
    theme_cowplot() +
    labs(x = NULL) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
}

create_collisson_boxplot <- function(metadata){
  md <- fix_data(metadata)
  ggplot(data = md, mapping= aes(x = collisson_clusters,
                                 y = ABSOLUTE.Purity)) +
    geom_boxplot(aes(colour = collisson_clusters), staplewidth = 0.5,
                 outliers = FALSE
    ) +
    geom_point(mapping = aes(x = collisson_clusters, y = ABSOLUTE.Purity),
               position = position_jitter(width = 0.1, height = 0.1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    labs(x = NULL)
  
}

create_tumor_fraction_plot <- function(metadata){
  md <- fix_data(metadata)
  # Rename columns and pivot the dataframe
  lookup <- c("Pathology_Review" = "Pathologist.Reviewed.Tumor.Cellularity",
              "Initial_slide_cellularity" = "Initial.Slide.Tumor.Cellularity",
              "Purity_class" = "Purity.Class..high.or.low.",
              "Leucocyte" = "DNA.methylation.leukocyte.percent.estimate",
              "ABSOLUTE" = "ABSOLUTE.Purity")
  
  md |> select(c(Pathologist.Reviewed.Tumor.Cellularity, Initial.Slide.Tumor.Cellularity, ABSOLUTE.Purity,Purity.Class..high.or.low.)) |>
    rename(any_of(lookup)) |>
    mutate(ABSOLUTE = ABSOLUTE*100) |>
    pivot_longer(cols = c(Pathology_Review,ABSOLUTE),
                 names_to = "method",
                 values_to = "cellularity"
    ) |>
    ggplot() +
    geom_boxplot(aes(x = method, y = cellularity),
                 staplewidth = 0.5,
                 outliers = FALSE) +
    geom_point(aes(y = cellularity, x = method, fill = Purity_class, color = Purity_class), position = position_jitterdodge()) +
    theme_cowplot() +
    labs(x = NULL, y = "Estimated Tumor Percentage")+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
  
  
}


create_radial_plot <- function(metadata){
  md <- fix_data(metadata)
  # Rename columns
  lookup <- c("Pathology_Review" = "Pathologist.Reviewed.Tumor.Cellularity",
              "Initial_slide_cellularity" = "Initial.Slide.Tumor.Cellularity",
              "Purity_class" = "Purity.Class..high.or.low.",
              "Leucocyte" = "DNA.methylation.leukocyte.percent.estimate")
  data <- md |>
    select(c(sample_id,bailey_clusters,moffit_clusters,collisson_clusters,Purity.Class..high.or.low., DNA.methylation.leukocyte.percent.estimate)) |>
    rename(any_of(lookup)) |>
    mutate(Purity_class = as.character(Purity_class)) |>
    arrange(bailey_clusters, collisson_clusters, moffit_clusters)
  # Define the number of samples
  n <- dim(data)[1]  
  data$sample <- paste0("S", 1:n)
  
  # Define colors for each category
  # Make sure to define in the same order as factor levels otherwise legends will be messed up later
  bailey_colors <- c("ADEX" = "magenta", "Immunogenic" = "red",
                     "Progenitor" = "royalblue", "Squamous" = "goldenrod" )
  
  collisson_colors <- c("Exocrine-like" = "green","Classical" = "cyan",
                        "Quasi-mesenchymal" = "deeppink")
  
  moffitt_colors <- c("Basal-like" = "orange", "Classical" = "darkblue")
  
  purity_colors <- c( "low" = "honeydew", "high" = "black")
  
  # Define spacing for groups
  # Define spacing: small gaps between samples, large gaps between groups
  group_counts <- table(data$bailey_clusters)  # Number of samples per group
  gaps <- rep(1, n)  # Default small gap for each sample
  
  # Add a larger gap after the last sample in each group
  cumulative_counts <- cumsum(as.numeric(group_counts))
  gaps[cumulative_counts] <- 6  # Increase space after groups
  # Add a gap at the start and end
  gaps[length(gaps)] <- 50  # Large gap at the end
  
  ## Set graphics parameters so we can display the figure and legend at the end
  #grid.newpage()
  #circle_size = unit(0.7, "npc") # snpc unit gives you a square region#
  
  #pushViewport(viewport(x = 0.4, y = 0.5, width = circle_size , height = circle_size))
  
  par(mfrow = c(1, 2))
  
  # Clear any existing plots
  circos.clear()
  
  # Set circular plot parameters with adjusted gaps
  circos.par(start.degree = 330, gap.after = gaps, cell.padding = c(0, 0, 0, 0),
             clock.wise = TRUE)
  
  # Initialize circular layout
  circos.initialize(factors = data$sample, xlim = c(0, 1))
  
  # Function to add annotation tracks
  # Function to add annotation tracks (without borders)
  add_circos_track <- function(data_column, color_mapping, track_height) {
    circos.track(
      factors = data$sample, ylim = c(0, 1),
      panel.fun = function(x, y) {
        sample_name <- get.cell.meta.data("sector.index")
        value <- data[data$sample == sample_name, data_column]
        circos.rect(0, 0, 1, 1, col = color_mapping[value], border =color_mapping[value] )  # Remove border
      },
      track.height = track_height
    )
}

# Add classification tracks
add_circos_track("bailey_clusters", bailey_colors, 0.1)
add_circos_track("collisson_clusters", collisson_colors, 0.1)
add_circos_track("moffit_clusters", moffitt_colors, 0.1)

# Add purity as an inner track
add_circos_track("Purity_class", purity_colors, 0.05)

# Add leukocyte score as a continuous gradient track
circos.track(
  factors = data$sample, ylim = c(0, 1),
  panel.fun = function(x, y) {
    sample_name <- get.cell.meta.data("sector.index")
    value <- data$Leucocyte[data$sample == sample_name]
    #print(value)
    col_pallette <- colorRampPalette(c("white", "red"))(100)[as.integer(value * 99 + 1)]
    circos.rect(0, 0, 1, 1, col = col_pallette, border = col_pallette)
  },
  track.height = 0.05
)

# Add legend to the plot
lgd_bailey = Legend(at = unique(data$bailey_clusters), title = "Bailey", legend_gp = gpar(fill = bailey_colors))
lgd_collisson = Legend(at = unique(data$collisson_clusters), title = "Collisson", legend_gp = gpar(fill = collisson_colors))
lgd_moffitt = Legend(at = unique(data$moffit_clusters), title = "Moffitt", legend_gp = gpar(fill = moffitt_colors))
lgd_purity = Legend(at = unique(data$Purity_class), title = "Purity(ABSOLUTE)", legend_gp = gpar(fill = purity_colors))
circos.clear()
col_fun <- colorRamp2(seq(0,100,length.out = 100),colorRampPalette(c("white", "red"))(100) )
lgd_leukocyte = Legend(col_fun = col_fun, title = "Leukocyte", at= c(0, 50, 100),
                       labels = c("0%", "50%", "100%"), direction = "horizontal")
pd = packLegend(lgd_purity, lgd_leukocyte,lgd_moffitt, lgd_collisson, lgd_bailey)
#draw(pd, x = unit(0.2, "npc") , y = unit(0.8, "npc"), 
#     just = c("top", "right"))
draw(pd)
par(mfrow = c(1, 1))
}
#create_radial_plot2(metadata)
#create_radial_plot(metadata)
#create_bailey_purity_boxplot(metadata)
#create_collisson_boxplot(metadata)
#create_moffitt_boxplot(metadata)
#create_tumor_fraction_plot(metadata)
