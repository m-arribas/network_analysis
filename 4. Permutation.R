# Load packages (make sure these are installed before loading)----------------------------------------------

library(tidyverse)
library(qgraph)
library(psychonetrics)
library(mice)
library(imputeTS)
library(viridis)
library(reshape2)


# Run the 3 models with the shared nodes and then output matrix variances from each network
# Follow the same steps (cleaning, development and visualisations) using the three seperate samples of interest 
# In our case, these are called PSY (psychotic disorders), bmd (bipolar disorders) and UMD (unipolar mood disorders) 
# These three models have the same nodes

network_data_umd <- read.csv("residuals_umd.csv")
network_data_bmd <- read.csv("residuals_bmd.csv")
network_data_psy <- read.csv("residuals_psy.csv")

# Define model parameters (these are the same across the 3 models)
network_data <- network_data_umd[, grepl("aggression|agitation|anxiety|cannabis|cognitive_impairment|disturbed_sleep|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|irritability|mood|paranoia|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking", names(network_data_umd) ) ]
net_var <- colnames(network_data) 

n_t        <- 6 # time points
n_v        <- 16 # number of symptom variables 
design_mat <- matrix(net_var, nrow = n_v, ncol = n_t, byrow = T) # design matrix for network & regression
rownames(design_mat) <- substr(design_mat[,1], 1, nchar(design_mat[, 1]) - 2)

# 1. Run the original 3 models ------------------------------------------------------------------------------

model_psy <- panelgvar(
  data = network_data_psy, # data
  vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements.
  estimator = 'ML',
  storedata = T
)

model_bmd <- panelgvar(
  data = network_data_bmd, # data
  vars = design_mat, 
  estimator = 'ML',
  storedata = T
)

model_umd <- panelgvar(
  data = network_data_umd, # data
  vars = design_mat, 
  estimator = 'ML',
  storedata = T
)

# Set seed for reproducible results
set.seed(1234)

# Run models
model_features_psy <- model_psy %>% 
  runmodel() %>% 
  fit()
model_features_bmd <- model_bmd %>% 
  runmodel() %>% 
  fit()
model_features_umd <- model_umd %>% 
  runmodel() %>% 
  fit()

# 2. Extract the adjacency matrix for temporal network ------------------------------------------------------------------------------

features_temporal_psy <- getmatrix(model_features_psy, "beta")
features_temporal_bmd <- getmatrix(model_features_bmd, "beta")
features_temporal_umd <- getmatrix(model_features_umd, "beta")

features_list <- list(features_temporal_psy, features_temporal_bmd, features_temporal_umd)

labelsAbb <- list("AGGR", "AGIT", "ANX", "CANN", "COGN", 
               "TOB", "SLEEP", "HALL", "INS", "IRR", 
               "MOOD", "PAR", "CONC", "MOTIV", "SUIC", 
               "TEAR")

for (p in 1:length(features_list)) {
  filename <- paste0(obs_list_strings[[p]], ".csv")
  data <- features_list[[p]]
  rownames(data) <- labelsAbb
  colnames(data) <- labelsAbb
  # sort columns and rows alphabetically
  data <- data %>% round(3)
  sorted_col_names <- sort(colnames(data))                                
  sorted_row_names <- sort(rownames(data)) 
  sorted_matrix <- data[sorted_row_names,sorted_col_names]                            
  write.csv(sorted_matrix, filename)
}

# 3. Plot the three temporal networks ------------------------------------------------------------------------------

g1_psy <- qgraph(features_temporal_psy,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )
g1_bmd <- qgraph(features_temporal_bmd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )
g1_umd <- qgraph(features_temporal_umd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )

# 4. Plot the centrality plots ------------------------------------------------------------------------------

qgraph::centrality(g1_psy)
qgraph::centrality(g1_bmd)
qgraph::centrality(g1_umd)

# 5. Calculate the observed difference in temporal edge estimates  ------------------------------------------------------------------------------

obs_diff_temp_umd_bmd <- features_temporal_umd - features_temporal_bmd
obs_diff_temp_umd_psy <- features_temporal_umd - features_temporal_psy
obs_diff_temp_bmd_psy <- features_temporal_bmd - features_temporal_psy

# save(obs_diff_temp_umd_bmd, obs_diff_temp_umd_psy, obs_diff_temp_bmd_psy, file = "/Users/Arribas/Desktop/PhD/Ongoing\ projects/Network/NA\ Code/Git/final/obs_diff_matrices.RData")

# 6. Permutation test ------------------------------------------------------------------------------

perm_diff_temp_umd_bmd <- list()
perm_diff_temp_umd_psy <- list()
perm_diff_temp_bmd_psy <- list()

cent_diff_temp_umd_bmd <- list()
cent_diff_temp_umd_psy <- list()
cent_diff_temp_bmd_psy <- list()

N_perm <- 250 
list_net_impute_scaled <- list()
start_time <- Sys.time()

network_data <- read.csv("total_2yr.csv") # Synthetic Data in Git

set.seed(1234)
network_data$index_diagnostic_group <- sample(c("UMD", "BMD", "PSY"), nrow(network_data), replace = TRUE)

for (i in 1:N_perm) {  # This will take a long time
  
    perm_data <- network_data %>%
    # Select columns to be sampled
    select(matches("aggression|agitation|anxiety|cannabis|cognitive_impairment|disturbed_sleep|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|irritability|mood|paranoia|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking")) %>%
    # Sample data
    mutate(across(everything(), ~ sample(.))) %>%
    # Bind with non-sampled demographic data
    bind_cols(network_data %>%
                select(matches("brcid_a|index_diagnostic_group|AgeF0|GenderF0|EthnicityF0|BoroughF0|ADF0|MSF0|ANXF0|ASYF0|count|total_words"))) 
      
    # Create subsets for different diagnostic groups
      perm_data_psy <- perm_data %>% filter(index_diagnostic_group == "PSY") %>% select(-index_diagnostic_group)
      perm_data_bmd <- perm_data %>% filter(index_diagnostic_group == "BMD") %>% select(-index_diagnostic_group)
      perm_data_umd <- perm_data %>% filter(index_diagnostic_group == "UMD") %>% select(-index_diagnostic_group)
  
  list_dfs <- c("perm_data_psy", "perm_data_bmd","perm_data_umd")  # make this a list we can subset through
  
  for (p in 1:3) {
    temporary_perm_data <- get(list_dfs[p])  
    
    # Linear imputation for numerical variables (not categorical variables: ethnicity, gender)
    
    impute_linear <- temporary_perm_data %>%
      select(-c('brcid_a', 'AgeF0', 'GenderF0', 'EthnicityF0','BoroughF0','ADF0', 'MSF0', 'ANXF0', 'ASYF0')) %>%
      na_interpolation() %>%
      round(0) %>%
      cbind(temporary_perm_data[, "brcid_a"]) %>% 
      setNames(c(colnames(.)[1:108], "brcid_a"))
    
    # Non-linear imputation for categorical variables only (ethnicity, gender)
    
    impute_non_linear <- temporary_perm_data %>%
      .[, c('AgeF0', 'GenderF0', 'EthnicityF0','BoroughF0','ADF0', 'MSF0', 'ANXF0', 'ASYF0')] %>%
      mutate(
        BoroughF0 = as.factor(BoroughF0),
        GenderF0 = as.factor(GenderF0),
        EthnicityF0 = as.factor(EthnicityF0),
        ADF0 = as.factor(ADF0),
        MSF0 = as.factor(MSF0),
        ANXF0 = as.factor(ANXF0),
        ASYF0 = as.factor(ASYF0)
      )
    
    impute_non_linear <- mice(impute_non_linear, method = c("", "pmm", "pmm", "", "", "", "", "")) %>%
      complete(., 2) %>%
      cbind(temporary_perm_data[, c("brcid_a")]) %>%
      setNames(c(colnames(.)[1:8], "brcid_a"))
    
    # Merge results for linear and nonlinear imputation
    
    net_impute <- merge(impute_linear, impute_non_linear, by = 'brcid_a') %>% as.data.frame()
    
    ## Regress out age, gender, medication vars, total words and score (number of EHR entries in FUP)
    
    for (k in 1:n_v){
      for (m in 1:n_t){
        x <- names(net_impute[design_mat[k, m]])
        y1<- names(net_impute[paste('total_words.', m, sep = '')])
        y2<- names(net_impute[paste('count.', m, sep = '')])
        f <- paste(x, "~", paste('AgeF0', 'GenderF0', 'EthnicityF0', 'ADF0', 'MSF0', 'ANXF0', 'ASYF0',y1 ,y2 , sep = ' + '))
        resid <- lm(f, net_impute) 
        net_impute[, x] <- resid$residuals
      }
    }
    
    # Center, scale & detrend
    
    net_impute <- net_impute %>% 
      mutate(across(net_var, ~ scale(.x, center = T, scale = T)))
    
    list_net_impute_scaled[[p]] <- net_impute
  }
  
  # Estimate the permuted networks
  
  model_psy <- panelgvar(
    data = list_net_impute_scaled[[1]], # data
    vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements. 
    estimator = 'ML',
    storedata = T
  )
  
  model_bmd <- panelgvar(
    data = list_net_impute_scaled[[2]], # data
    vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements. 
    estimator = 'ML',
    storedata = T
  )
  
  model_umd <- panelgvar(
    data = list_net_impute_scaled[[3]], # data
    vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements. 
    storedata = T
  )
  
  #Run models
  perm_model_features_psy  <- model_psy %>% 
    runmodel()
  perm_model_features_bmd  <- model_bmd %>% 
    runmodel()
  perm_model_features_umd  <- model_umd %>%
    runmodel()
  
  # Extract the adjacency matrix - temporal
  perm_features_temporal_psy <- getmatrix(perm_model_features_psy, "beta") %>% 
    as.data.frame() %>%
    round(3)
  perm_features_temporal_bmd <- getmatrix(perm_model_features_b,d, "beta") %>% 
    as.data.frame() %>%
    round(3)
  perm_features_temporal_umd <- getmatrix(perm_model_features_umd, "beta") %>% 
    as.data.frame() %>%
    round(3)
  
  perm_diff_temp_umd_bmd[[i]] <- perm_features_temporal_umd - perm_features_temporal_bmd
  perm_diff_temp_umd_psy[[i]] <- perm_features_temporal_umd - perm_features_temporal_psy
  perm_diff_temp_bmd_psy[[i]] <- perm_features_temporal_bmd - perm_features_temporal_psy
  
  # Add centrality

  labelsAbb <- list("AGGR", "AGIT", "ANX", "CANN", "COGN", 
                    "TOB", "SLEEP", "HALL", "INS", "IRR",
                    "MOOD", "PAR", "CONC", "MOTIV", "SUIC", 
                    "TEAR")
  
  g1_psy <- qgraph(perm_features_temporal_psy,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )
  g1_bmd <- qgraph(perm_features_temporal_bmd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )
  g1_umd <- qgraph(perm_features_temporal_umd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.0 )
   
  cent_psy <-qgraph::centrality(g1_psy)
  cent_bmd <-qgraph::centrality(g1_bmd)
  cent_umd <-qgraph::centrality(g1_umd)
  
  cent_perm_diff_umd_bmd_temp <- cent_umd - cent_bmd
  cent_perm_diff_umd_psy_temp <- cent_umd - cent_psy
  cent_perm_diff_bmd_psy_temp <- cent_bmd - cent_psy
  
}

# save(perm_diff_temp_umd_bmd, perm_diff_temp_umd_psy, perm_diff_temp_bmd_psy, file = "/Users/Arribas/Desktop/PhD/Ongoing\ projects/Network/NA\ Code/Git/final/perm_diff_matrices.RData")

#7.  Calculate the corrected p-values ------------------------------------------------------------------------------

names_list_temp <- list("temp_umd_bmd", "temp_umd_psy", "temp_bmd_psy")

obs_list_temp <- list(obs_diff_temp_umd_bmd, obs_diff_temp_umd_psy, obs_diff_temp_bmd_psy)

perm_list_temp <- list(list_perm_diff_temp_umd_bmd, list_perm_diff_temp_umd_psy, list_perm_diff_temp_bmd_psy)

# Loop over each edge in the matrix

for (p in 1:length(obs_list_temp)) {
  filename <- paste0(names_list_temp[p], ".csv")
  filename_2 <- paste0(names_list_temp[p], "_corrected.csv")
  
  # Initialise a matrix to store p-values
  
  p_values_matrix <- matrix(NA, nrow = 16, ncol = 16)
  rownames(p_values_matrix) <- labelsAbb
  colnames(p_values_matrix) <- labelsAbb
  
  # Number of permutations
  n_perm <- length(perm_list_temp[[p]])
  
  for (i in 1:16) {
    for (j in 1:16) {
      
      # Extract observed diff for this edge
      obs_diff <- obs_list_temp[[p]][i, j]
      
      # Calculate how many permuted differences are as extreme or more extreme
      extreme_count <- sum(sapply(perm_list_temp[[p]], function(x) abs(x[i, j]) >= abs(obs_diff)))  
      
      # Calculate the p-value and store it in the matrix 
      p_values_matrix[i, j] <- extreme_count / n_perm %>%
        round(3)
      
      # sort columns and rows alphabetically
      sorted_col_names <- sort(colnames(p_values_matrix))                                
      sorted_row_names <- sort(rownames(p_values_matrix)) 
      
      p_values_matrix <- p_values_matrix[sorted_row_names, sorted_col_names]                            
      
      # write.csv(p_values_matrix, filename)   
      
      # Correction for multiple comparisons                            
      flat_p_values <- as.vector(p_values_matrix)
      adjusted_p_values <- p.adjust(flat_p_values, method= "fdr") %>%
        round(3)
      adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow=16, ncol=16)
      
      rownames(adjusted_p_values_matrix) <- labelsAbb
      colnames(adjusted_p_values_matrix) <- labelsAbb
      
      # sort columns and rows alphabetically
      sorted_col_names <- sort(colnames(adjusted_p_values_matrix))                                
      sorted_row_names <- sort(rownames(adjusted_p_values_matrix)) 
      
      adjusted_p_values_matrix <- adjusted_p_values_matrix[sorted_row_names,sorted_col_names]                            
      
      write.csv(adjusted_p_values_matrix, filename_2)
      
    }
  }
}


# 7.  Make a directional matrix ------------------------------------------------------------------------------

#  Lets make a direction matrix for the effect sizes, and only display the results with a corrected p-value <.05

corrected_p_vals_list <- c("temp_umd_bmd_corrected.csv","temp_umd_psy_corrected.csv","temp_bmd_psy_corrected.csv")

for (p in 1:length(obs_list_temp)) {
  filename <- paste0(names_list_temp[p], "_sig_directions.csv")
  corrected_p_vals <- read.csv(corrected_p_vals_list[p], row.names = 1) %>%
    as.matrix()
  
  # Calculate the mean of the permuted edge weights
  permuted_means <- apply(simplify2array(perm_list_temp[[p]]), 1:2, mean)
  # Initialise a matrix to store the direction of differences
  direction_matrix <- matrix(NA, nrow = 16, ncol = 16)
  
  for (i in 1:16) {
    for (j in 1:16) {
      # Extract observed diff for this edge
      obs_diff <- obs_list_temp[[p]][i, j]
      
      if(obs_diff > permuted_means[i, j]) {
        direction_matrix[i, j] <- "positive"
      } else if(obs_diff < permuted_means[i,j]) { 
        direction_matrix[i, j] <- "negative"
      } else {
        direction_matrix[i, j] <- "neutral" 
      }
    }
    
  }
  rownames(direction_matrix) <- labelsAbb
  colnames(direction_matrix) <- labelsAbb
  
  # sort columns and rows alphabetically
  sorted_col_names <- sort(colnames(direction_matrix))                                
  sorted_row_names <- sort(rownames(direction_matrix)) 
  direction_matrix <- direction_matrix[sorted_row_names, sorted_col_names]     
  
  # Only display the direction that have a significant p-val
  significance_matrix <- corrected_p_vals < 0.05
  numeric_direction_matrix <- matrix(0, nrow = 16, ncol = 16)
  numeric_direction_matrix[direction_matrix == "positive"] <- 1
  numeric_direction_matrix[direction_matrix == "negative"] <- -1
  numeric_direction_matrix[direction_matrix == "neutral"] <- 0
  result_matrix <- significance_matrix * numeric_direction_matrix
  
  write.csv(result_matrix, filename)
  
}

# 7.  Plot histograms only for the comparisons that had signficant results------------------------------------------------------------------------------

# Given the actual data ---------------------------------------------------------------

load("corrected_p_vals_matrices.RData")     # Corrected p-values - temp_umd_bmd_corrected, temp_umd_psy_corrected, temp_bmd_psy_corrected

load("sig_directions_matrices.RData")   # Significant directions - temp_umd_bmd_sig_directions, temp_umd_psy_sig_directions, temp_bmd_psy_sig_directions

load("temp_obs_diff_matrices.RData")   # Observed differences - obs_diff_temp_umd_bmd, obs_diff_temp_umd_psy, obs_diff_temp_bmd_psy

plot_heatmap_temp <- function(comparison) { 
  if (comparison == "temp-umd-psy") {
    title <- "Temporal heatmap for UMD-PSY"
    filename <- paste0("/path/to/save/file/",comparison, ".png")  # Change this to the directory where you want to save images
    effect_sizes <-  obs_diff_temp_umd_psy
    directions <- temp_umd_psy_sig_directions
    p_values <- temp_umd_psy_corrected }
  
  else if (comparison == "temp-umd-bmd") {
    title <- "Temporal heatmap for UMD-BMD"
    filename <- paste0("/path/to/save/file/",comparison, ".png")
    effect_sizes <-  obs_diff_temp_umd_bmd
    directions <- temp_umd_bmd_sig_directions
      p_values <- temp_umd_bmd_corrected}
  
    else if (comparison == "temp-bmd-psy") {
    title <- "Temporal heatmap for BMD-PSY"
    filename <- paste0("/path/to/save/file/",comparison, ".png")
    effect_sizes <-  obs_diff_temp_bmd_psy
    directions <- temp_bmd_psy_sig_directions
      p_values <- temp_bmd_psy_corrected }
    
  rownames(effect_sizes) <- labelsAbb
  colnames(effect_sizes) <- labelsAbb
  
  # sort columns and rows alphabetically
  sorted_col_names <- sort(colnames(effect_sizes))                                
  sorted_row_names <- sort(rownames(effect_sizes)) 
  effect_sizes <- effect_sizes[sorted_row_names, sorted_col_names] %>%
    reshape2::melt()  
  
  directions <-  directions %>% 
    reshape2::melt()  
  
  p_values <- p_values %>% 
    reshape2::melt()  
  
  combined_data <- cbind (effect_sizes, Direction = directions$value, P_value = p_values$value)
  colnames(combined_data)[1:3] <- c("From", "To", "Effect")
  combined_data$Sig <- combined_data$P_value < 0.05
  
  # Create heatmap
  ggplot(combined_data, aes(x = To, y = From, fill = Effect)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = ifelse(Sig, "*", "")),
              size = 8, vjust = 0.7, colour = "white") +
    scale_fill_viridis_c(limit = c(-0.075, 0.075), 
                         breaks = c(-0.075, -0.050, -0.025, 0, 0.025, 0.050, 0.075), 
                         labels = c("-0.075", "-0.050", "-0.025", "0","0.025", "0.050", "0.075"),
                         name = "Effect Size") +
    theme_minimal() + 
    coord_fixed() +
    theme(axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(size = 15) ) +
    xlab("to") +
    ylab("from") + 
    ggtitle(title)
  ggsave(filename, dpi = 600)
}


# Comparison UMD-PSY

plot_heatmap_temp("temp-umd-psy") 
plot_heatmap_temp("temp-umd-bmd") 
plot_heatmap_temp("temp-bmd-psy") 

