# Load packages (make sure these are installed before loading)----------------------------------------------

library(tidyverse)
library(gtsummary)
library(data.table)
library(mice)
library(imputeTS)
library(caret)
library(standardize)

# Feature variables ------------------------------------------------------------------------------

# Load data as CSV
total_cleaned <- read.csv("/path/to/file/dynamic_syn.csv") # Synthetic Data in Git

cols_symptoms <- c(
  'aggression', 'agitation', 'anergia', 'anhedonia', 'anxiety', 
  'apathy', 'appetite', 'arousal', 'auditory_hallucination', 
  'bad_dreams', 'blunted_flat_affect', 'circumstantial_speech', 
  'cognitive_impairment', 'concrete_thinking', 'delusion', 
  'derailment_of_speech', 'disturbed_sleep', 'diurnal',
  'early_morning', 'echolalia', 'elation', 'emotional_withdrawn',
  'flight_of_ideas', 'formal_thought_disorder', 'grandiosity', 
  'guilt', 'hallucination', 'hallucination_otg', 'helpless', 
  'hopeless', 'hostility', 'insight', 'insomnia', 'irritability',
  'lonely', 'loss_of_coherence', 'low_energy', 'mood_instability', 
  'mutism', 'negative_symptom', 'nightmare', 'paranoia', 'passivity',
  'persecutory', 'poor_concentration', 'poor_motivation', 
  'poverty_of_speech', 'poverty_of_thought', 'social_withdrawal', 
  'stupor', 'suicidal', 'tangential_speech', 'tearful', 
  'thought_block', 'thought_broadcast', 'thought_insert', 
  'thought_withdrawal', 'visual_hallucination', 'waxy_flexibility',
  'weightloss', 'worthless', 'cannabis', 'cocaine', 'mdma', 
  'current_smoking')

cols_additional_variables <- c(
  'brcid_a', "words", 'index_diagnostic_group', 'month_index')

total_2yr <- total_cleaned %>%
  select(all_of(cols_additional_variables), all_of(cols_symptoms))

# Filter the data to obtain EHR entries that are 2 years prior to the index date

total_2yr <- total_2yr %>%
  filter(month_index >= -24 & month_index < -6 )

# Generate Follow Up Point "fup" variable (2-year follow up period is split into six 3mo windows between T-24mo and T-6mo which we call FUPs)

total_2yr <- total_2yr %>%
  mutate(fup = case_when(
    month_index >= -24 & month_index < -21 ~ 1,
    month_index >= -21 & month_index < -18 ~ 2,
    month_index >= -18 & month_index < -15 ~ 3,
    month_index >= -15 & month_index < -12 ~ 4,
    month_index >= -12 & month_index < -9  ~ 5,
    month_index >= -9  & month_index < -6  ~ 6,
    TRUE ~ NA_real_
  ))

total_2yr_fu <- total_2yr %>%
  group_by(brcid_a, fup) %>%
  summarise(count = n(), .groups= 'drop', total_words = sum(words)) # The new var "count" denotes the number of EHR entries per FUP

# Make all the symptom variables binary (0/1) within each 3mo time window

total_2yr <- total_2yr %>% 
  group_by(brcid_a, fup) %>%
  summarise(across(cols_symptoms, sum)) %>% 
  data.table() %>% # Convert to data.table for higher efficiency with large datasets
  mutate(across(all_of(cols_symptoms), ~ as.integer(. > 0)))

# Calculate mean and SD number of follow-up points

nfups <- total_2yr %>%
  group_by(brcid_a) %>%
  summarise(count=n())

mean(nfups$count)
sd(nfups$count)

# Remove patients with 4 or less FUPs

exclusions <- nfups %>% 
  filter(count == 1 | count == 2 | count == 3 | count == 4) # these are the IDs we will exclude

# Keep individuals with more than 4 FUPs

total_2yr <- total_2yr %>% 
  filter(!(brcid_a %in% exclusions$brcid_a))  %>% 
  inner_join(total_2yr_fu, by=c("brcid_a", "fup")) %>% 
  reshape(idvar = "brcid_a", timevar = "fup", direction = "wide")  

# Order all the column names

# Separate 'brcid_a' from the other column names
other_colnames <- setdiff(colnames(total_2yr), "brcid_a")
ordered_other_colnames <- other_colnames[order(other_colnames)]
ordered_colnames <- c("brcid_a", ordered_other_colnames)

total_2yr <- total_2yr[, ..ordered_colnames]

# Demographic/Clinical variables ------------------------------------------------------------------------------

demog <- total_cleaned %>%
  select(c("brcid_a", "index_diagnostic_group", "age_index", "Gender_ID", 
           "ethnicitycleaned", "borough", 'ad_new_index', 'ms_new_index', 
           'anx_new_index', 'asy_new_index')) %>% 
  filter(brcid_a %in% total_2yr$brcid_a) %>%
  unique()

total_2yr <- # Load the filtered population from the data cleaning above
  inner_join(demog, total_2yr, by = "brcid_a") %>% 
  as.data.frame() %>%
  unique() #  Merge to the demographics

# Create Demographics summary table

demog_tbl_summary <-
  total_2yr %>%
  mutate(ethnicitycleaned = factor(ethnicitycleaned, levels = c(
    "White", "Black", "Other", "Asian", "Mixed", "missing"))) %>%
  select(age_index, Gender_ID, ethnicitycleaned, borough, index_diagnostic_group) %>%
  tbl_summary(
    label = list(age_index ~ "Age", 
                 Gender_ID ~ "Gender", 
                 ethnicitycleaned ~ "Ethnicity",
                 borough ~ "Borough"),
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"),
    missing = "no" ,
    digits = age_index ~ c(1, 1)
  ) %>%
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "demographics_table.tex")

# Create Clinical Variable summary table

clinical_tbl_summary <-
  total_2yr %>%
  select(ad_new_index, ms_new_index, anx_new_index, asy_new_index, index_diagnostic_group) %>%
  tbl_summary(
    type = list( "ad_new_index" ~ "categorical", 
                 "ms_new_index" ~ "categorical",
                 "anx_new_index" ~ "categorical", 
                 "asy_new_index" ~ "categorical"),
    label = list( ad_new_index ~ "Antidepressants",
                  ms_new_index ~ "Mood stabilisers",
                  anx_new_index ~ "Anxiolytics",
                  asy_new_index ~ "Antipsychotics"),
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"),
    missing = "ifany" ,
    digits = all_continuous() ~ 2, 
  ) %>%
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "clinical.tex")

# Create final dataset ------------------------------------------------------------------------------

# Rename the sociodemographic and clinical variables

total_2yr <- total_2yr %>%
  mutate(Gender_ID = recode(Gender_ID, 
                            'Female' = '2', 
                            'Male' = '1', 
                            'Other' = '0', 
                            "missing" = 'NA')) %>%
  mutate(ethnicitycleaned = recode(ethnicitycleaned, 
                                   'Asian' = '4',
                                   'Black' = '3',
                                   'Mixed' = '2',
                                   'White' = '1',
                                   'Other' = '0', 
                                   "missing" = 'NA')) %>%
  select(-c("index_diagnostic_group")) %>%
  rename(AgeF0 = age_index,
         GenderF0 = Gender_ID,
         EthnicityF0 = ethnicitycleaned,
         BoroughF0 = borough,
         ADF0 = ad_new_index,
         MSF0 = ms_new_index,
         ANXF0 = anx_new_index,
         ASYF0 = asy_new_index) %>% 
  mutate(across(c('GenderF0', 'EthnicityF0', "BoroughF0"), factor))

write.csv(total_2yr, 'total_2yr.csv') #  This will be re-used for R script "3. Bootstrapping.R"

# Linear imputation for numerical variables (not categorical variables: ethnicity, gender)

set.seed(1234)
impute_linear <- total_2yr %>%
  select(-c('brcid_a', 'AgeF0', 'GenderF0', 'EthnicityF0','BoroughF0','ADF0', 'MSF0', 'ANXF0', 'ASYF0')) %>%
  na_interpolation() %>%
  round(0) %>%
  cbind(total_2yr[, "brcid_a"]) %>%
  setNames(c(colnames(.)[1:402], "brcid_a"))

# Non-linear imputation for categorical variables only (ethnicity, gender)

impute_non_linear <- total_2yr %>%
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

impute_non_linear <- mice(impute_non_linear, method = c("", "pmm", "pmm", "", "", "", "", ""),  seed = 22112 ) %>%
  complete(., 2) %>%
  cbind(total_2yr[, c("brcid_a")]) %>%
  setNames(c(colnames(.)[1:8], "brcid_a"))

# Merge results for linear and nonlinear imputation

net_impute <- merge(impute_linear, impute_non_linear, by = 'brcid_a') %>% as.data.frame()

# 2. Remove Near-Zero Variance (NZV) variables -------------------------------------------------------------------------------

nzv_all <- net_impute %>% 
  select(-c("count.1", "count.2", "count.3", "count.4", "count.5", "count.6", 
            "total_words.1", "total_words.2", "total_words.3", "total_words.4", "total_words.5", "total_words.6",
            "BoroughF0", 'AgeF0', 'GenderF0', 'EthnicityF0', 'ADF0', 'MSF0', 'ANXF0', 'ASYF0')) %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  nearZeroVar(. , saveMetrics=TRUE)

# write.csv(nzv, "nzv_all.csv")

final_vars <- nzv_all %>%
  filter(nzv == FALSE) %>%
  rownames(.) %>% 
  unique()


# 3. Regression -------------------------------------------------------------------------------

# Create a logical vector using grepl
net_impute_vars <- net_impute[, grepl("aggression|agitation|anxiety|cannabis|cocaine|cognitive_impairment|delusion|disturbed_sleep|emotional_withdrawn|hopeless|guilt|\\bhallucination.1\\b|\\bhallucination.2\\b|\\bhallucination.3\\b|\\bhallucination.4\\b|\\bhallucination.5\\b|\\bhallucination.6\\b|hostility|irritability|mood|paranoia|poor_concentration|insight|poor_motivation|suicidal|tearful|current_smoking|weightloss", names(net_impute))]

# Set up design matrix

net_var <- colnames(net_impute_vars)
n_t        <- 6 # number of time points
n_v        <- 23 # number of symptom variables (after removing demographics variables)
design_mat <- matrix(net_var, nrow = n_v, ncol = n_t, byrow = T) # design matrix for network & regression
rownames(design_mat) <- substr(design_mat[, 1], 1, nchar(design_mat[, 1]) - 2)
design_mat

# Regress out age, gender, medication vars, total words and score (number of EHR entries in FUP) 

set.seed(1234)

for (k in 1:n_v){
  for (i in 1:n_t){
    x <- names(net_impute[design_mat[k,i]])
    y1 <- names(net_impute[paste('total_words.', i, sep = '')])
    y2 <- names(net_impute[paste('count.', i, sep = '')])
    f <- paste(x, "~", paste('AgeF0', 'GenderF0', 'EthnicityF0', 'ADF0', 'MSF0', 'ANXF0', 'ASYF0', y1, y2, sep = ' + '))
    resid <- lm(f, net_impute)
    net_impute[, x] <- resid$residuals
  }
}

# Center, scale & detrend

set.seed(1234)
net_impute <- net_impute %>% 
  mutate(across(net_var, ~ scale(.x, center = T, scale = T)))
dim(net_impute)

write.csv(net_impute, "residuals.csv") # This will be used in R script "2. Network.R"

# Sanity checks -

#  plotting detrended residuals (as aggregate rather than individual data)

ggplot(net_impute %>%
         pivot_longer(aggression.1 : aggression.6,  # Replace "aggression" for the name of the feature to plot
                      names_to = 'Time', values_to = 'Value'), 
       aes(Time, Value)) +
  geom_violin(trim = FALSE) +
  stat_summary(color = 'red', size = 0.2, fun.data = mean_cl_normal) +
  labs(title = 'Aggression', x = 'Time (FUP)') +
  scale_x_discrete(labels = c('FU1', 'FU2', 'FU3', 'FU4', 'FU5', 'FU6')) +
  theme_bw() +
  theme(text = element_text(size = 9))

