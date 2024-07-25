# Load necessary libraries
library(shellpipes)
library(randomForestSRC)
library(survival)
library(satpred)
library(readxl) # For reading Excel files
library(dplyr) # For data manipulation
satpredtheme()

set.seed(8888)

# Read data from Excel file
data_path <- "C:/Users/luyu/Desktop/ISLR/data_surf.xlsx"
train_df <- read_excel(data_path, skip = 2)

# List of binary variables
binary_vars <- c("ECTR", "HP", "HD", "CRRT", "PE", "HP_only", "sex", 
                 "mental_health_illness", "dementia", "asthma", 
                 "respiratory_failure", "respiratory_disease", 
                 "hypertension", "coronary_heart_disease", 
                 "congestive_heart_failure", "cardiac_disease", 
                 "peptic_ulcer", "abnormal_liver_function", 
                 "abnormal_kidney_function", "diabetes", "malignant_tumor", 
                 "sequelae_of_stroke", "rheumatic_immune_diseases", 
                 "history_of_past_mental_health_disease", "airway",  
                 "CPR", "ROSC", "intubation_or_ventilation", "intubation", 
                 "ventilation", "induced_vomiting_or_catharsis", 
                 "induced_vomiting", "catharsis", "special_antidote", 
                 "vasoactive_drugs", "glucocorticoids", "activated_carbon", 
                 "gastric_lavage", "gastric_lavage_special_bed", 
                 "gastric_lavage_comp", 
                 "endotracheal_intubation_during_gastric_lavage", 
                 "gastric_lavage_after_admission", 
                 "local_hospital_has_gastric_lavage_before_admission")

# List of numeric variables
numeric_vars <- c("dosage_liquid_1", "age", "MAP", "HR", "RR")

# Keep only necessary variables and remove missing values
train_df <- train_df %>%
  select(all_of(binary_vars), all_of(numeric_vars), 
         "survival_time", "censor") %>%
  na.omit()

# Convert necessary columns to factors
train_df <- train_df %>%
  mutate(across(all_of(binary_vars), as.factor))

# Remove columns with only one level and print their names
cols_with_one_level <- sapply(train_df, function(col) is.factor(col) && length(unique(col)) == 1)
cols_with_one_level_names <- names(cols_with_one_level[cols_with_one_level])
train_df <- train_df %>% select(-one_of(cols_with_one_level_names))

# Print column names with only one level
print(cols_with_one_level_names)

# Set up parameter grid for tuning
params_rfsrc <- expand.grid(mtry = c(2, 3), nodesize = seq(2, 20, length.out = 5))

# Fit Random Survival Forest (RSF) model
tuned_rfsrc <- modtune(Surv(survival_time, censor) ~ ., train_df, param_grid = params_rfsrc,
                       modfun = rfsrc.satpred, nodedepth = NULL, forest = TRUE, parallelize = TRUE, seed = 8888)

print(tuned_rfsrc$besTune)
plot(tuned_rfsrc)

# Fit the model with the best parameters
fit_rfsrc <- modfit(tuned_rfsrc, return_data = FALSE)

# Generate individual survival curves
scurves_rfsrc <- get_indivsurv(fit_rfsrc, train_df)
plot(scurves_rfsrc)

# Calculate concordance score
concord_rfsrc <- get_survconcord(fit_rfsrc)
print(concord_rfsrc)

# Compute permutation variable importance
vimp_rfsrc <- get_varimp(fit_rfsrc, type = "perm", newdata = train_df, nrep = 20, modelname = "rfsrc")
plot(vimp_rfsrc)

# Save results
saveVars(fit_rfsrc, scurves_rfsrc, concord_rfsrc, vimp_rfsrc)