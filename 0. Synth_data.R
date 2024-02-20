library(dplyr)
library(charlatan)

# Define the number of unique IDs
num_ids <- 10000

# Consistent columns
consistent_columns <- data.frame(
  brcid_a = 1:num_ids,
  index_diagnostic_group = sample(c("UMD", "BMD", "PSY"), num_ids, replace = TRUE),
  age_index = sample(18:70, num_ids, replace = TRUE),
  Gender_ID = sample(c("Male", "Female", "missing"), num_ids, replace = TRUE),
  ethnicitycleaned = sample(c("White", "Black", "Asian", "missing", "Other"), num_ids, replace = TRUE),
  borough = sample(c("A", "B", "C", "D"), num_ids, replace = TRUE),
  ad_new_index = sample(c("1", "0"), num_ids, replace = TRUE),
  ms_new_index =  sample(c("1", "0"), num_ids, replace = TRUE),
  anx_new_index =  sample(c("1", "0"), num_ids, replace = TRUE),
  asy_new_index =  sample(c("1", "0"), num_ids, replace = TRUE)
)


# Generate a random number of rows for each ID (between 5 and 20)
replication_factors <- sample(5:20, num_ids, replace = TRUE)

# Replicate consistent data
replicated_consistent_columns <- do.call(rbind, 
                                         Map(function(df, n) df[rep(1, n), ], 
                                             split(consistent_columns, consistent_columns$brcid_a), 
                                             replication_factors))
rownames(replicated_consistent_columns) <- NULL # Reset row names

varying_column_names <-c( 'aggression','agitation','anergia','anhedonia','anxiety',
                          'apathy','appetite','arousal','auditory_hallucination',
                          'bad_dreams','blunted_flat_affect','circumstantial_speech',
                          'cognitive_impairment','concrete_thinking','delusion',
                          'derailment_of_speech','disturbed_sleep','diurnal',
                          'early_morning','echolalia','elation','emotional_withdrawn',
                          'flight_of_ideas','formal_thought_disorder','grandiosity',
                          'guilt','hallucination','hallucination_otg','helpless',
                          'hopeless','hostility','insight','insomnia','irritability',
                          'lonely','loss_of_coherence','low_energy','mood_instability',
                          'mutism','negative_symptom','nightmare','paranoia','passivity',
                          'persecutory','poor_concentration','poor_motivation',
                          'poverty_of_speech','poverty_of_thought','social_withdrawal',
                          'stupor','suicidal','tangential_speech','tearful','thought_block',
                          'thought_broadcast','thought_insert','thought_withdrawal',
                          'visual_hallucination','waxy_flexibility','weightloss','worthless',
                          'cannabis','cocaine','mdma','current_smoking')

generate_numeric_column <- function(n) {
  sample(0:10, n, replace = TRUE)
}

# Function to generate all varying data
generate_varying_data <- function(n) {
  # Generate numeric columns
  numeric_columns <- setNames(
    lapply(varying_column_names, function(name) generate_numeric_column(n)),
    varying_column_names
  )
  
  # Generate categorical columns (replace 'Category1', 'Category2' with actual categories)
  words = sample(1:50000, n, replace = TRUE)
  month_index = sample(-30:6, n, replace = TRUE)
  
  # Combine all columns
  data_frame <- cbind(do.call(data.frame, numeric_columns), 
                      words = words, 
                      month_index = month_index)
  
  return(data_frame)
}

varying_columns <- generate_varying_data(nrow(replicated_consistent_columns))

# Combine consistent and varying data
synthetic_df <- cbind(replicated_consistent_columns, varying_columns)

# View the head of the dataframe
head(synthetic_df)

# Write to CSV
write.csv(synthetic_df, "dynamic_syn.csv")
