library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
#' 
#' read.delim("my_data.csv", sep = ",")
#' 
read_data <- function(intensity_data, delimiter) {
  data <- read.delim(intensity_data, sep = delimiter)
  #trans_data <- t(data)
  #scaled_trans_data <- scale(trans_data, center = TRUE, scale = TRUE)
  #scaled_data <- t(scaled_trans_data)
  return(as.data.frame(data))
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  sum_pca <- summary(pca_results)
  
  imp_sum_pca <- sum_pca$importance
  
  var_pca <- imp_sum_pca[2,]
  
  
  return(var_pca)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  sum_pca <- summary(pca_results)

  imp_sum_pca <- sum_pca$importance

  cml_var <- imp_sum_pca[3,] #3
  names_pca <- names(cml_var) #1
  tibbl <- tibble(variance_explained=pca_ve, 
                  principal_components=as.factor(names_pca), 
                  cumulative=cml_var)
  
  
  return(tibbl)
}

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble generated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples
plot_pca_variance <- function(variance_tibble) {
  pca_var_plot <- ggplot(data = variance_tibble, aes(x=reorder(principal_components,-variance_explained), y=variance_explained)) +
    geom_bar(stat='identity', color="black", aes(fill="Variance Explained"))+
    geom_line(aes(y=cumulative, group=1, colour="Cumulative"))+
    geom_point(aes(y=cumulative,  colour="Cumulative"))+
    labs(title = 'PCA Variance plot',
         x= 'PC', y='% Variance')+
    scale_colour_manual("Cumulative", breaks ="Cumulative",
                        values ="black") +
    scale_fill_manual("Variance Explained", breaks = "Variance Explained", values="lightblue")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.line = element_line(),
          panel.background = element_rect(fill = 'white'),
          plot.caption = element_text(hjust = 0),
          plot.caption.position = "plot",
          legend.box.just = c("left", "top"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          legend.key.size = unit(5, units = "mm"))
  
   return(pca_var_plot)
}
# Load dataset from github

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
  pca_x <- pca_results$x
  meta <- read.csv(metadata) %>%
    select(geo_accession, SixSubtypesClassification)
  pca_x_tbl <- as_tibble(pca_x, rownames='geo_accession')
  joined_tbl <- pca_x_tbl %>%
    left_join(meta, by='geo_accession')
  bi_plot <- ggplot(joined_tbl) + 
    geom_point(aes(x=PC1, y=PC2, color=SixSubtypesClassification)) +
    theme_classic() + 
    labs(x="PC1", y="PC2")
  return(bi_plot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_csv (str): The path to the differential expression results
#'   file we have provided
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the CSV.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
  diff_exp_data <- read.csv(diff_exp_csv)
  diff_exp_data_tibbl <- as_tibble(diff_exp_data, rownames = "probe_id")
  index_that_pass_thres <- which(diff_exp_data_tibbl$padj < fdr_threshold)
  sig_probes <- diff_exp_data_tibbl$probe_id[index_that_pass_thres]
  return(sig_probes)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  intensity_tibbl <- as_tibble(intensity, rownames="probe_id")
  filetered_intensity_tibbl <- filter(intensity_tibbl, probe_id %in% sig_ids_list)
  matrix_data <- as.matrix(column_to_rownames(filetered_intensity_tibbl, var='probe_id'))
  return(matrix_data)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  col_pal <- colorRampPalette(brewer.pal(num_colors, "RdBu"))(num_colors)
  heatmap_plot <- heatmap(de_intensity, col=col_pal)
  return(heatmap_plot)
}