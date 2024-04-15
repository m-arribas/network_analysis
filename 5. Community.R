rm(list=ls())

# Libraries all installed now just need to load ---------------------------------------------------------------

library(tidyverse)
library(bootnet)
library(psychonetrics)
library(qgraph)
library(mgm)
library(graphicalVAR)
library(tidyverse)
library(mlVAR)
library(foreign)
library(mice)
library(Hmisc)
library(imputeTS)
library(patchwork)
library(ggridges)
library(ordinal)
library(standardize)
library(ggsci)
library(igraph)
library(optimr)
library(reshape2)

# Community analysis detection 

# Load the data  ---------------------------------------------------------------

# Run the 3 models with the shared nodes and then output matrix variances from each network
# In our case, these are called PSY (psychotic disorders), bmd (bipolar disorders) and UMD (unipolar mood disorders) 
# These three models have the same nodes 
# To run the community analysis on the SMD model, follow the same steps but using residuals from the SMD model and adapting the number of nodes from 16 to 23. The spin for the community analysis should also be modified from 8 to 12.

network_data_umd <- read.csv("residuals_umd.csv")  # Found in the data folder
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

labelsAbb <- list("AGGR", "AGIT", "ANX", "CANN", "COGN", 
                  "TOB", "SLEEP", "HALL", "INS", "IRR", 
                  "MOOD", "PAR", "CONC", "MOTIV", "SUIC", 
                  "TEAR")

# 3. Extract the three temporal networks as igraphs ------------------------------------------------------------------------------

g1_psy <- qgraph(features_temporal_psy,  labels = labelsAbb,  diag = TRUE,  threshold = 0.01 )
g1_bd <- qgraph(features_temporal_bd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.01 )
g1_nbmd <- qgraph(features_temporal_nbmd,  labels = labelsAbb,  diag = TRUE,  threshold = 0.01 )

g_psy <- as.igraph(g1_psy, attributes = TRUE)
g_bd <- as.igraph(g1_bd, attributes = TRUE)
g_nbmd <- as.igraph(g1_nbmd, attributes = TRUE)

# 4. Create 1000 community structures with Spinglass algorithm ------------------------------------------------------------------------------

model <- g_psy #  OR g_bd OR g_nbmd 
spin <- vector("list", 1000)

for (i in 1:1000) {
     set.seed(i)
     spin[[i]] <- spinglass.community(model, spin=8, implementation = "neg") 
     }

# 5. Create a covariance matrix to visualise results ------------------------------------------------------------------------------

n=16
co_occ_matrix <- matrix(0, 16, 16)

for (communities in spin) {
    membership <- communities$membership
    for (i in 1:n) {
        for (j in 1:n) {
            if (membership[i] == membership[j]) {
                co_occ_matrix[i,j]  <- co_occ_matrix[i,j] +1
                }
            }
        }
    }

co_occ_matrix <- co_occ_matrix / 1000
# write.csv(co_occ_matrix, "covariance_matrix.csv") 
# saveRDS(spin, "community_list_1000_it_restricted_spin8.RDS")


# 6. Visualise covariance matrix as a heatmap ------------------------------------------------------------------------------

#  write.csv(top_100, "top_100_communities_restricted_spin8_neg.csv")

co_occ_matrix <- co_occ_matrix[-1]
row.names(co_occ_matrix) <- labelsAbb
colnames(co_occ_matrix) <- labelsAbb
co_occ_matrix <- as.matrix(co_occ_matrix)

hc <- hclust(as.dist(1-co_occ_matrix))
plot(hc, hang = -1)
dendro <- as.dendrogram(hc)
order <- order.dendrogram(dendro)

reordered_co_occ_matrix <- co_occ_matrix[, order]

reordered_co_occ_matrix <- reordered_co_occ_matrix[order ,]

reordered_co_occ_matrix[reordered_co_occ_matrix == 1] <- 0

long_co_occ_matrix <- melt(reordered_co_occ_matrix)

long_co_occ_matrix$value <- as.numeric(long_co_occ_matrix$value)

x <- ggplot(long_co_occ_matrix, aes(Var1, Var2, fill= value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_distiller(name = "Covariance strength",palette = "Blues", direction = +1, limits = c(0,1)) +
  theme_minimal() +
  xlab("Node1")+
  ylab("Node2")+
  theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),axis.text.x = element_text(angle=90, hjust=1, size=15),axis.text.y = element_text(size=15) ) +
  ggtitle("Co-occurence heatmap")

# ggsave("covariance_matrix_reordered.pdf", x)


# 7. Find the top most occurring individual communities ------------------------------------------------------------------------------

# Function to get a standardized string representation of each community
community_to_string <- function(community) {
  split_communities <- split(1:length(community$membership), community$membership)
  sapply(split_communities, function(x) paste(sort(x), collapse = "-"))
}

# Assuming 'spin' is your list of 1000 communities objects, convert each community to a standardized string format
community_strings <- lapply(spin, community_to_string)

# Flatten the list and count occurrences
community_tabulation <- table(unlist(community_strings))

# Sort and get the names of top 3 most frequent communities
top_communities <- sort(community_tabulation, decreasing = TRUE)
top_communities <- names(sort(community_tabulation, decreasing = TRUE))
top_100 <- top_communities[1:100]

#  write.csv(top_100, "top_100_communities_restricted_spin8_neg.csv")