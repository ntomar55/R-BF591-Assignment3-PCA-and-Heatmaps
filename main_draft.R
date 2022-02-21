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
  trans_data <- t(data)
  scaled_trans_data <- scale(trans_data, center = TRUE, scale = TRUE)
  scaled_data <- t(scaled_trans_data)
  return(as.data.frame(scaled_data))
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



res <- read_data("example_intensity_data.csv", " ")
fake_pca_results1 <- prcomp(res, center = FALSE, scale = FALSE)
fake_pca_results1
sum_pca <- summary(fake_pca_results1)

imp_sum_pca <- sum_pca$importance

var_pca <- imp_sum_pca[2,]
var_pca


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
  tibbl <- tibble(names_pca, pca_ve, cml_var)
  
  
  return(tibbl)
}


imp_sum_pca
colnames(imp_sum_pca)
rownames(imp_sum_pca)
var_pca

colnames(var_pca)
names(var_pca)
rownames(var_pca)
var_pca
cml_var <- imp_sum_pca[3,] #3
names_pca <- names(cml_var) #1

tbl <- tibble(names_pca, var_pca, cml_var)
tbl

data_frame1 <- as.data.frame(tbl)


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
  pca_var_plot <- ggplot(data = variance_tibble, aes(x=reorder(names_pca,-var_pca), y=var_pca)) +
    geom_bar(stat='identity', color="black", aes(fill="Variance Explained"))+
    geom_line(aes(y=cml_var, group=1, colour="Cumulative"))+
    geom_point(aes(y=cml_var,  colour="Cumulative"))+
    labs(title = 'PCA Variance plot',
         caption = "The bar chart represents individual variances of each principal component.
The connected dots represent the cumulative variances",
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



tbl
tbl$names_pca
data_frame1$names_pca

data_frame1$var_pca
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
  pca_rot <- scale(pca_results$rotation,  center=T, scale=T)
  meta <- read.csv(metadata) %>%
    select(geo_accession, SixSubtypesClassification)
  pca_rot_tbl <- as_tibble(pca_rot_data, rownames='geo_accession')
  joined_tbl <- pca_rot_tbl %>%
    left_join(metadata, by='geo_accession')
  bi_plot <- ggplot(joined_tbl) + 
    geom_point(aes(x=-PC1*100, y=PC2*100, color=SixSubtypesClassification)) +
    theme_classic() + 
    labs(x="PC1", y="PC2")
  return(bi_plot)
}

as_tibble(fake_pca_results1$rotation, rownames='aa')

dim(fake_pca_results1$x)
imp_sum_pca
pca_res_rot <- as.data.frame(scale(fake_pca_results1$rotation[,1:2], center=TRUE, scale=TRUE))


ggplot(data = pca_res_rot, aes(x = PC1, y = PC2)) +
  geom_point()+
  labs(title = 'PCA plot', x= "PC1", y="PC2")

test_metadata <- read.csv("test.csv")
test_metadata
metadata <- read.csv('proj_metadata.csv') %>%
  select(geo_accession, SixSubtypesClassification)
names(metadata)
dim(metadata)

pca_rot_data <- scale(fake_pca_results1$rotation, center=T, scale=T)
#pca_rot_data <- fake_pca_results1$rotation
tbl2 <- as_tibble(pca_rot_data, rownames='geo_accession')
tbl2
joined_data <- tbl2 %>%
  left_join(metadata, by='geo_accession')

joined_data$SixSubtypesClassification

ggplot(joined_data) + 
  geom_point(aes(x=-PC1*100, y=PC2*100, color=SixSubtypesClassification)) +
  theme_classic() + 
  labs(x="PC1", y="PC2")

metadata$geo_accession

fake_pca_results1$x %>% 
  as_tibble(rownames='geo_accession')
as_tibble(pca_data, rownames='aaa')
labeled <- fake_pca_results1$x %>% 
  as_tibble(rownames='geo_accession') %>% 
  left_join(metadata, by='geo_accession')
labeled$SixSubtypesClassification

labeled %>% 
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, color=SixSubtypesClassification)) +
  theme_classic()

dim(labeled)
head(labeled)  
rownames(tbl2)
metadata$geo_accession
dim(metadata)
dim(res)


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

test_diff_exp <- read.csv('test.csv')
tbl3 <- as_tibble(test_diff_exp, rownames='probe_id')
#idx <- which(tbl3$padj < 0.01)
#idx
#sig <- tbl3$probe_id[idx]
#sig
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
fx_intensity_df <- data.frame(GSM971958 = runif(1000, 0, 15),
                              GSM971959 = runif(1000, 0, 15),
                              GSM971967 = runif(1000, 0, 15),
                              GSM971968 = runif(1000, 0, 15))
row.names(fx_intensity_df) <- paste0(1:1000, "_s_at")
fx_intensity_mat <- as.matrix(fx_intensity_df)

test_mat <- fx_intensity_mat[c('1_s_at', '2_s_at'),]
function_mat <- return_de_intensity(fx_intensity_df, c('1_s_at', '2_s_at'))
function_mat

tibl4 <- as_tibble(fx_intensity_df, rownames="probe_id")
id_lists <- list_significant_probes('test.csv', 0.01)
id_lists
f_tibl4 <- filter(tibl4, probe_id %in% id_lists)
f_tibl4 
a <- as.matrix(column_to_rownames(f_tibl4, var='probe_id') )
a
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
  
  return()
}

sig_p_ids <- list_significant_probes('differential_expression_results.csv', 0.01)
sig_p_ids
inten_mat <- return_de_intensity(res, sig_p_ids)
inten_mat
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(11)
heatmap(inten_mat, col=coul, scale="column")