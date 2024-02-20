# Load packages (make sure these are installed before loading)----------------------------------------------

library(tidyverse)
library(psychonetrics)
library(qgraph)
library(mgm)
library(graphicalVAR)
library(mlVAR)
library(foreign)
library(mice)
library(Hmisc)
library(imputeTS)
library(ggridges)
library(ordinal)
library(standardize)
library(ggsci)

# Load data

network_data <- read.csv("/Users/Arribas/Desktop/PhD/Ongoing\ projects/Network/NA\ Code/Git/final/data/residuals.csv")  # Output from 1. Pre-processing.R script 

network_data <-  network_data[, grepl("aggression|agitation|anxiety|cannabis|cocaine|cognitive_impairment|delusion|disturbed_sleep|emotional_withdrawn|hopeless|guilt|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|hostility|irritability|mood|paranoia|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking|weightloss", names(network_data))]

# 1. Model development -------------------------------------------------------------------------------

n_t        <- 6 # time points
n_v        <- 23 # number of symptom variables (after removing demographics variables)
net_var <- colnames(network_data) 
design_mat <- matrix(net_var, nrow = n_v, ncol = n_t, byrow = T) # design matrix for network & regression
rownames(design_mat) <- substr(design_mat[,1], 1, nchar(design_mat[, 1]) - 2)
design_mat


model <- panelgvar(
  data = network_data, # network data
  vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements (follow up point)
  estimator = 'ML',
  storedata = F
)

# Run models

# Saturated model
set.seed(1234)
model_features_sat <- model %>%
  runmodel() %>%
  fit()

# Sparse model
set.seed(1234)
model_features_spa <- model_features_sat %>%
  prune(alpha = 0.01, adjust = "none") %>% 
  fit()

# Compare if saturated or sparse model performs better
psychonetrics::compare(
  saturate = model_features_sat,
  sparse   = model_features_spa
)

# saveRDS(model_features_sat,'model_features_sat.RData')

features_temporal <- getmatrix(model_features_sat, "beta")
features_contemporaneous <- getmatrix(model_features_sat, "omega_zeta_within")
features_between <- getmatrix(model_features_sat, "omega_zeta_between")

features_covariances <- list(temp = features_temporal, contemp = features_contemporaneous, between = features_between)

# saveRDS(features_covariances,'features_covariances.RData')

# 2. Estimating Model Recovery -------------------------------------------------------------------------------

# Generate data from network
set.seed(1234)
generate_data <- model_features_sat %>% 
  generate(7041) # Make this of same size as actual population size used for model 

model_rec <- panelgvar(
  data = generate_data, # data generated
  vars = design_mat, # The same design matrix, with a row indicating a variable and a column a wave of measurements.
  estimator = 'ML'
)

# Run recovery model
set.seed(1234)
model_features_rec <- model_rec %>% 
  runmodel() %>%
  fit()

#saveRDS(model_features_rec,'model_features_rec.RData')

# 3. Visualising networks ------------------------------------------------------------------------------

features_covariances <- readRDS("/Users/Arribas/Desktop/PhD/Ongoing\ projects/Network/NA\ Code/Git/final/data/features_covariances.RData")

# Extract networks:
features_temporal <- features_covariances[[1]]
features_contemporaneous <- features_covariances[[2]]
features_between <- features_covariances[[3]]

labels     <- rownames(design_mat) %>% as.list()

labelsAbb <- list("AGGR", "AGIT", "ANX", "CANN", "COC", "COGN", 
                  "TOB", "DEL", "SLEEP", "EMOT", "GUIL", "HALL", 
                  "HOPE", "HOST", "INS", "IRR", "MOOD", "PAR", 
                  "CONC", "MOTIV", "SUIC", "TEAR", "WGHT")

# pdf("network_temporal.pdf", width =100, height=100)
g1 <- qgraph(features_temporal, layout = "spring", repulsion = 0.82, mar = c(5, 5, 5, 5),
             labels = labelsAbb, curve = 1, curveAll = TRUE, 
             diag = TRUE, title = "Temporal", 
             vsize = 8, label.prop = 1, label.scale = TRUE, label.scale.equal = TRUE, # vsize = size of nodes, label.cex = scalar on label size 
             asize = 3, threshold = 0.022,   # mar = margin, asize = size of arrowhead, thresholc = removed from network,
             edge.label.cex = 1.2, edge.label.position = 0.45, esize=7, #  edge.label.cex = size of edge weight
             theme = "colorblind") 
box("figure")
dev.off()

# pdf("network_contempo.pdf", width =100, height=100)
g2 <- qgraph(features_contemporaneous, layout = g1$layout, mar=c(5, 5, 5, 5),
             labels = labelsAbb, curve = 0.6, curveAll = TRUE, 
             diag = FALSE, title = "Contemporaneous",  
             vsize = 8, label.prop = 1, label.scale = TRUE, label.scale.equal= TRUE, 
             asize = 3, threshold = 0.08,   
             edge.label.cex = 1.2, edge.label.position = 0.5, esize = 7, 
             theme = "colorblind" ) #
box("figure")
dev.off()

# pdf("network_between.pdf", width =100, height=100)
g3 <- qgraph(features_between, layout = g1$layout, mar=c(5,5,5,5), 
             labels = labelsAbb, curve = -0.6, curveAll = TRUE,  
             diag = FALSE, title = "Between-individuals",
             vsize = 8, label.prop = 1, label.scale= TRUE,label.scale.equal= TRUE,
             asize = 3, threshold = 0.17,  
             edge.label.cex = 1.2, edge.label.position = 0.55, esize = 7, 
             theme = "colorblind") #
box("figure")
dev.off()


# 4. Centrality Figures  ------------------------------------------------------------------------------


x <-qgraph::centrality(g1)
in_temp <- x[1] %>% 
  as.data.frame() %>%
  round(3)

out_temp <- x[2] %>% 
  as.data.frame() %>% 
  round(3)

centr_contemp <- qgraph::centrality(g2)[1]  %>%
  as.data.frame() %>% 
  round(3)

centr_betw <- qgraph::centrality(g3)[1]  %>%
  as.data.frame() %>%
  round(3)


total <- in_temp %>%
  merge(out_temp, by = 'row.names', all = TRUE) %>%
  rename(
    Feature = Row.names,
  ) %>%
  merge(centr_contemp, by.x = "Feature", by.y = 'row.names', all = TRUE) %>%
  rename(  OutDegree = OutDegree.x, Contemporaneous = OutDegree.y ) %>%
  merge(centr_betw, by.x = "Feature", by.y = 'row.names', all = TRUE) %>%
  rename(  OutDegree = OutDegree.x, Between = OutDegree.y ) 


# write.csv(total,"network_SMD_centrality.csv")

# First, sort the dataframe by centrality
total <- total[order(total$OutDegree), ]

# Create an additional column for OutDegree the nodes
total$order <- seq_along(total$OutDegree)

# Now create the ggplot (can repeat this for temporal_in, contemprananeous and between-subejct)

# png("network_centrality_temp_out.png", res=150, height=1000, width=400)
ggplot(total, aes(x = order, y = OutDegree)) +
  geom_point() +  # Add points
  geom_line() +   # Connect the points with lines
  scale_x_continuous(breaks = total$order, labels = total$Feature) +  # Adjust x-axis to show node labels
  scale_y_continuous(breaks = seq(0,0.13,0.02)) +  # Adjust y-axis limits
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_blank()) +
  ylab("centrality (out)") + 
  coord_flip() # Rotate x-axis labels for readability

#dev.off()