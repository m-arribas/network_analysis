# Load packages (make sure these are installed before loading)----------------------------------------------

library(tidyverse)
library(ggplot2)
library(psychonetrics)
library(qgraph)
library(mgm)
library(graphicalVAR)
library(mlVAR)
library(foreign)
library(mice)
library(imputeTS)
library(ggridges)
library(ordinal)
library(standardize)
library(ggsci)
library(bootnet)


# 1. Bootstrapping loop -------------------------------------------------------------------------------

# Load the data 

total_2yr <- read.csv("total_2yr.csv") # outputted from "1.Pre-processing.R" script

network_data <- total_2yr[, grepl("aggression|agitation|anxiety|cannabis|cocaine|cognitive_impairment|delusion|disturbed_sleep|emotional_withdrawn|hopeless|guilt|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|hostility|irritability|mood|paranoia|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking|weightloss", names(total_2yr) )]
dim(network_data) 
net_var <- colnames(network_data) 
n_t        <- 6 # time points
n_v        <- 23 # number of symptom variables (after removing demographics variables)
design_mat <- matrix(net_var, nrow = n_v, ncol = n_t, byrow = T) # design matrix for network & regression
rownames(design_mat) <- substr(design_mat[, 1], 1, nchar(design_mat[, 1]) - 2)
design_mat

# Create empty lists to store results from the bootstrapping loop

features_boot_temporal <- list()
features_boot_contempo <- list()
features_boot_betweens <- list()
av_fitness             <- list() # average goodness of fit

##This will take a long time

start_time <- Sys.time()
reps <- 250 # Number of repeats

for (i in 1:reps) {
  defaultW <- getOption("warn")
  options(warn = -1) # ignores warning sign
  rowstouse        <- base::sort(sample(nrow(total_2yr), nrow(total_2yr)*.75)) #This takes 75% of pop
  ni_data          <- total_2yr[rowstouse,]
  
  ## Linear interpolation of all numeric variables (except sociodemog and meds)
  
  net <- select(ni_data, -c('X', 'brcid_a', 'AgeF0', 'GenderF0', 'EthnicityF0', 
                            'BoroughF0', 'ADF0', 'MSF0', 'ANXF0', 'ASYF0'))
  net_impute <- round(na_interpolation(net), 0)   
  net_impute_linear <- cbind(net_impute, ni_data[, c(2)])
  colnames(net_impute_linear)[403] <- "brcid_a" # checked
  
  ## Non-linear interpolation of  gender + ethnicity
  
  nonlinear_matrix <- ni_data[, c('AgeF0', 'GenderF0', 'EthnicityF0', 'BoroughF0', 
                                  'ADF0', 'MSF0', 'ANXF0', 'ASYF0')] %>%
    mutate(
      BoroughF0 = as.factor(BoroughF0),
      GenderF0 = as.factor(GenderF0),
      EthnicityF0 = as.factor(EthnicityF0)
    )
  
  net_impute_nonlinear <- mice(nonlinear_matrix, method = c("", "pmm", "pmm", "", "", "", "", ""))
  
  completedata <- complete(net_impute_nonlinear, 2) %>%
    select(GenderF0, EthnicityF0) %>%
    cbind(ni_data[, c("brcid_a")]) 
  colnames(completedata)[3] <- "brcid_a"
  
  net_impute <- merge(net_impute_linear, completedata, by = 'brcid_a')
  
  cv_data <- ni_data %>%
    select(brcid_a, AgeF0, ADF0, MSF0, ANXF0, ASYF0) %>%
    merge(net_impute, by = 'brcid_a') %>%
    drop_na() %>%
    mutate(
      GenderF0 = as.factor(GenderF0),
      EthnicityF0 = as.factor(EthnicityF0),
      ADF0 = as.factor(ADF0),
      MSF0 = as.factor(MSF0),
      ANXF0 = as.factor(ANXF0),
      ASYF0 = as.factor(ASYF0)
    )
  
  ## Regress out age, gender, medication vars, total words and score (number of EHR entries in FUP)
  for (k in 1:n_v){
    for (t in 1:n_t){
      x <- names(cv_data[design_mat[k, t]])
      y1<- names(cv_data[paste('total_words.', t, sep = '')])
      y2<- names(cv_data[paste('count.', t, sep = '')])
      f <- paste(x, "~", paste('AgeF0', 'GenderF0', 'EthnicityF0', 'ADF0', 'MSF0', 'ANXF0','ASYF0', y1, y2, sep = ' + '))
      resid <- lm(f, cv_data)
      cv_data[, x] <- resid$residuals
    }
  }
  
  # Center, scale & detrend
  net_impute_scaled <- cv_data %>%
    mutate(across(net_var, ~ scale(.x, center = T, scale = T)))
  net_impute_scaled <- net_impute_scaled[, grepl("aggression|agitation|anxiety|cognitive_impairment|delusion|disturbed_sleep|emotional_withdrawn|hopeless|guilt|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|hostility|irritability|mood|paranoia|persecutory|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking|weightloss", names(net_impute_scaled))]
  
  # Set up bootstrapped model
  model_boot <- panelgvar(
    data = net_impute_scaled, # data
    vars = design_mat, # The design matrix, with a row indicating a variable and a column a wave of measurements. Note that NA indicates missing variables
    estimator = 'ML', 
    storedata = F
  )
  
  # Run bootstrapped model
  features_boot <- model_boot %>%
    runmodel()
  features_boot_temporal[[i]] <- getmatrix(features_boot, "beta")
  features_boot_contempo[[i]] <- getmatrix(features_boot, "omega_zeta_within")
  features_boot_betweens[[i]] <- getmatrix(features_boot, "omega_zeta_between")
  av_fitness[[i]] <- features_boot %>% 
    fit() %>% 
    as.data.frame()
  
  print(cat(paste('\n Completed ', i, ' of ', reps)))
  options(warn = defaultW)
  
}
end_time <- Sys.time()
total_time <- end_time - start_time
total_time

#saveRDS(features_boot_temporal, "features_boot_temporal.RDS")
#saveRDS(features_boot_contempo, "features_boot_contempo.RDS")
#saveRDS(features_boot_betweens, "features_boot_betweens.RDS")
#saveRDS(av_fitness, "av_fitness.RDS")


# 2. Averaging post-bootsrap -------------------------------------------------------------------------------


# Define functions  ------------------------------------------------------------------------------

# 95% CIs
lCI    <- function(mean, sd, reps = reps){(mean - ((1.96 * (sd/sqrt(reps)))))}
uCI    <- function(mean, sd, reps = reps){(mean + ((1.96 * (sd/sqrt(reps)))))}

# Clean the list
sortMe <- function(x, n_v, reps = 250){
  avTempMat      <- do.call(cbind, x)
  avTempMatX     <- array(avTempMat, dim=c(dim(x[[1]]), length(x)))
  avTempMatm     <- apply(avTempMatX, c(1, 2), mean, na.rm = T)
  avTempMatsd    <- apply(avTempMatX, c(1, 2), sd, na.rm = T)
  avTempMatcil   <- avTempMatm
  avTempMatciu   <- avTempMatm
  avTempMatse    <- avTempMatm
  avTempMatm_sig <- avTempMatm
  for (i in 1:n_v){
    for (j in 1:n_v){
      avTempMatcil[i, j] <- lCI(avTempMatm[i, j], avTempMatsd[i, j], reps)
      avTempMatciu[i, j] <- uCI(avTempMatm[i, j], avTempMatsd[i, j], reps)
      if (avTempMatcil[i, j] < 0 & avTempMatciu[i, j] > 0){
        avTempMatm_sig[i, j] <- 0
      }
    }
  }
  
  avTempMatse <- avTempMatsd/sqrt(reps)
  
  return(
    list(
      mean = avTempMatm_sig,
      se   = avTempMatse,
      upperCI = avTempMatciu,
      lowerCI = avTempMatcil)
  )
}

# This is the only data we can provide

features_boot_temporal <- readRDS("features_boot_temporal.RDS") %>%
  unlist(recursive = FALSE)
features_boot_contempo <- readRDS("features_boot_contempo.RDS") %>%
  unlist(recursive = FALSE)
features_boot_betweens <- readRDS("features_boot_betweens.RDS") %>%
  unlist(recursive = FALSE)

av_temporal <- sortMe(features_boot_temporal, 23, 250)  # specify number of features (23) + number of iterations (250)
av_contempo <- sortMe(features_boot_contempo, 23, 250)
av_betweens <- sortMe(features_boot_betweens, 23, 250)

n_v = 23
n_t = 6
labelsD = list("AGGR", "AGIT", "ANX", "CANN", "COC", 
               "COGN", "TOB", "DEL", "SLEEP", "EMOT",
               "GUIL", "HALL", "HOPE", "HOST", "INS",
               "IRR", "MOOD", "PAR", "CONC", "MOTIV", 
               "SUIC", "TEAR", "WGHT")


# Adapt for the network of interest (example shown here for temporal network, but can adapt to contemporaneous/between-subject)

# Find upper confidence interval

upper <- av_temporal$upperCI %>% 
  as.data.frame() %>%
  rename('aggression' = 1, 'agitation' = 2, 'anxiety'= 3, 'cannabis use'= 4, 'cocaine use'= 5, 
         'cognitive impairment' = 6, 'tobacco use' = 7, 'delusional thinking' = 8, 'disturbed sleep' = 9, 'emotional withdrawal' = 10,
         'guilt' = 11, 'hallucinations (all)' = 12, 'feeling hopeless' = 13, 'hostility' = 14, 'poor insight' = 15, 
         'irritability' = 16, 'mood instability' = 17, 'paranoia' = 18, 'poor concentration' = 19, 'poor motivation' = 20, 
         'suicidality' = 21,'tearful' = 22, 'weightloss' = 23) %>% 
  mutate(V2 = labelsD)   %>% 
  pivot_longer(1:all_of(n_v), names_to = 'V1', values_to = 'upperCI')  

# Find lower confidence interval

lower <- av_temporal$lowerCI %>% 
  as.data.frame() %>%
  rename('aggression' = 1, 'agitation' = 2, 'anxiety'= 3, 'cannabis use'= 4, 'cocaine use'= 5, 
         'cognitive impairment' = 6, 'tobacco use' = 7, 'delusional thinking' = 8, 'disturbed sleep' = 9, 'emotional withdrawal' = 10,
         'guilt' = 11, 'hallucinations (all)' = 12, 'feeling hopeless' = 13, 'hostility' = 14, 'poor insight' = 15, 
         'irritability' = 16, 'mood instability' = 17, 'paranoia' = 18, 'poor concentration' = 19, 'poor motivation' = 20, 
         'suicidality' = 21,'tearful' = 22, 'weightloss' = 23) %>% 
  mutate(V2 = labelsD)   %>% 
  pivot_longer(1:all_of(n_v), names_to = 'V1', values_to = 'lowerCI')  

# Merge the lower and upper confidence intervals - if crosses zero - force to 0

allCI <- merge(upper, lower) %>%
  mutate(do_not_cross_zero = ifelse( ((upperCI > 0 & lowerCI < 0) |  (upperCI < 0 & lowerCI> 0)), 0, 1)) %>% 
  pivot_wider(names_from='V1', values_from = 'do_not_cross_zero') %>% 
  as.data.frame() %>% 
  select(-c('upperCI', 'lowerCI'))

allCI[is.na(allCI)] <- 0

allCI <- allCI %>%
  group_by(V2) %>% 
  summarise(across(aggression:weightloss, sum)) %>% 
  select(-c('V2'))


correlation_coeff_boot_means <- av_temporal$mean %>%
  as.data.frame() %>%
  round(3) %>%
  rename('aggression' = 1, 'agitation' = 2, 'anxiety'= 3, 'cannabis use'= 4, 'cocaine use'= 5, 
         'cognitive impairment' = 6, 'tobacco use' = 7, 'delusional thinking' = 8, 'disturbed sleep' = 9, 'emotional withdrawal' = 10,
         'guilt' = 11, 'hallucinations (all)' = 12, 'feeling hopeless' = 13, 'hostility' = 14, 'poor insight' = 15, 
         'irritability' = 16, 'mood instability' = 17, 'paranoia' = 18, 'poor concentration' = 19, 'poor motivation' = 20, 
         'suicidality' = 21,'tearful' = 22, 'weightloss' = 23) 

correlation_coeff_boot_SE <- av_temporal$se %>%
  as.data.frame() %>%
  round(4) %>%
  rename('aggression' = 1, 'agitation' = 2, 'anxiety'= 3, 'cannabis use'= 4, 'cocaine use'= 5, 
         'cognitive impairment' = 6, 'tobacco use' = 7, 'delusional thinking' = 8, 'disturbed sleep' = 9, 'emotional withdrawal' = 10,
         'guilt' = 11, 'hallucinations (all)' = 12, 'feeling hopeless' = 13, 'hostility' = 14, 'poor insight' = 15, 
         'irritability' = 16, 'mood instability' = 17, 'paranoia' = 18, 'poor concentration' = 19, 'poor motivation' = 20, 
         'suicidality' = 21,'tearful' = 22, 'weightloss' = 23) 

mean <- correlation_coeff_boot_means*allCI
SE <- correlation_coeff_boot_SE*allCI

mean[] <- paste0(as.matrix(mean), " (", as.matrix(SE), ")")
rownames(mean) <- c('aggression', 'agitation', 'anxiety', 'cannabis use', 'cocaine use', 
                    'cognitive impairment', 'tobacco use' ,'delusional thinking', 'disturbed sleep', 'emotional withdrawal', 
                    'guilt', 'hallucinations (all)', 'feeling hopeless', 'hostility', 'poor insight', 
                    'irritability', 'mood instability', 'paranoia', 'poor concentration', 'poor motivation', 
                    'suicidality', 'tearful', 'weightloss')

# write.csv(mean, "mean_SE_bootstrapped_av_temporal.csv")

# 3. Plotting bootstrapped vs actual model estimates  ------------------------------------------------------------------------------

# Load actual model estimates

features_covariances <- readRDS('features_covariances.RData') # This is the only data that we can provide (not individual level)
features_temporal <- features_covariances[[1]]
features_contemporaneous <- features_covariances[[2]]
features_between <- features_covariances[[3]]

intervals <- av_temporal$mean %>% 
  as.data.frame() %>%
  rename('AGGR' = 1, 'AGIT' = 2, 'ANX' = 3, 'CANN' = 4, 'COC' = 5,
         'COGN' = 6, 'TOB' = 7, 'DEL' = 8, 'SLEEP' = 9, 'EMOT' = 10,
         'GUIL' = 11, 'HALL' = 12, 'HOPE' = 13, 'HOST' = 14, 'INS' = 15,
         'IRR' = 16, 'MOOD' = 17, 'PAR' = 18, 'CONC' = 19, 'MOTIV' = 20,
         'SUIC' = 21, 'TEAR' = 22, 'WGHT' = 23) %>%
  mutate(V2 = labelsD)  %>% 
  pivot_longer(1:all_of(n_v), names_to='V1', values_to = 'Correlation')  %>%
  arrange(Correlation) %>%
  mutate(CI_U = av_temporal$upperCI %>%    # upper limit averaged from bootstrapping
           as.data.frame() %>%
           pivot_longer(1:all_of(n_v), names_to= 'V1', values_to = 'Correlation') %>%
           arrange(Correlation) %>%
           dplyr::select(Correlation) %>%
           unlist(),
         CI_L = av_temporal$lowerCI %>%    # lower limit averaged from bootstrapping
           as.data.frame() %>%
           pivot_longer(1:all_of(n_v), names_to= 'V1', values_to = 'Correlation') %>%
           arrange(Correlation) %>%
           dplyr::select(Correlation) %>%
           unlist(),
         RealCorr = features_temporal %>%  # whole data estimates
           as.data.frame() %>%
           pivot_longer(1:all_of(n_v), names_to ='V1', values_to = 'Correlation') %>%
           arrange(Correlation) %>%
           dplyr::select(Correlation) %>%
           unlist(),
         Same = ifelse(V1 == V2, '#DA7422', 'black'),
         AlphaCorr  = ifelse((CI_U > 0 & CI_L < 0) |
                               (CI_U < 0 & CI_L > 0), 0.001, 0.1),
         V1 = paste(V2, '-', V1))


intervals_order <- intervals[order(abs(intervals$Correlation), decreasing=TRUE), ] 
intervals_order$Correlation <- abs(intervals_order$Correlation)
intervals_order$RealCorr <- abs(intervals_order$RealCorr)
intervals_order$CI_U <- abs(intervals_order$CI_U)
intervals_order$CI_L <- abs(intervals_order$CI_L)

autocorrelations <- ggplot(intervals_order[,] ) +
  geom_point(aes(Correlation, reorder(V1, Correlation), alpha = AlphaCorr, color = 'Bootstrap')) +
  geom_pointrange(aes(Correlation, reorder(V1, Correlation), xmin = CI_U, xmax = CI_L, color = 'Bootstrap'), alpha = 0.5) +   #scale_x_discrete(labels=c("0.05"=".05","0.06"=".06","0.07"=".07","0.08"=".08","0.09"=".09","0.10"=".10"))+
  geom_point(aes(RealCorr, V1, color = 'Full Data'), alpha = 1) +
  scale_alpha(guide = 'none', range = c(0.05, 0.75)) +
  scale_color_manual(name = 'Network', values = c('black', '#D64933')) +
  labs(x = 'Estimate') +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title.x = element_text(size = 14))

 ggsave("/Users/Arribas/Desktop/PhD/Ongoing\ projects/Network/NA\ Code/Git/temporal_bootstrapped_vs_actual_estimates.png", width = 4)
