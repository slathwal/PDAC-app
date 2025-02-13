# Use data - rppa
# Also use the metadata

# Samples have to be in rows and variables in columns
# So subjects in rows and protein measurements in columns
library(factoextra)
library(cowplot)

create_pca_plot <- function(data, metadata, col_variable){
  
  metadata <- metadata |> mutate(Tumor.Sample.ID = str_replace_all(metadata$Tumor.Sample.ID, "-","."))
  
  #sample_metadata <- merge(samples, metadata, by.x = "sample_id", by.y = "Tumor.Sample.ID")
  
  res.pca <- prcomp(data, scale = FALSE)
  summary(res.pca)
  variance_explained <- res.pca$sdev^2/sum(res.pca$sdev^2)
  names(variance_explained) <- colnames(res.pca$x)
  
  # plot the scree plot using fviz_eig from factoextra package

  fviz_eig(res.pca, addlabels = TRUE) + 
    labs(title = "Scree Plot",
         x = "Principal component",
         y = "% of variances explained")
  
  samples <- data.frame(sample_id = rownames(data))
  samples <- data.frame(samples, res.pca$x)
  pca_sample_metadata <- merge(samples, metadata, by.x = "sample_id", by.y = "Tumor.Sample.ID")
  
  # Plot PC1 and PC2 with legends for samples (include metadata as well)
  
  pca_sample_metadata |>
    ggplot(aes(x=PC1, y = PC2)) +
    geom_point(aes(color = .data[[col_variable]]), size = 4) +
    labs(x=paste0("PC1: ",round(variance_explained["PC1"]*100,1),"%"),
         y=paste0("PC2: ",round(variance_explained["PC2"]*100,1),"%"))+
    theme_cowplot()
}


#pca_sample_metadata <- merge(samples, metadata, by.x = "sample_id", by.y = "Tumor.Sample.ID")



#proteomics_subset <- proteomics[rownames(proteomics) %in% sample_metadata[sample_metadata$Purity.Class..high.or.low. == "high",]$sample_id,]






#samples <- data.frame(sample_id = rownames(proteomics_subset))
#sample_metadata <- merge(samples, metadata, by.x = "sample_id", by.y = "Tumor.Sample.ID")


# create_pca_plot(proteomics, metadata, "gender")

