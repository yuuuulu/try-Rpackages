# method 1--using ranger package(quicker and not exploded R)
# Load necessary packages
library(ranger)
library(survival)
library(dplyr)
library(ggplot2)

# Load the lung dataset from the survival package
data("lung")

# Define binary variables and numeric variables
binary_vars <- c("sex", "ph.ecog")

numeric_vars <- c("age", "ph.karno", "pat.karno", "meal.cal", "wt.loss")

# Create a new dataframe containing only the specified variables
data1 <- lung %>%
  select(all_of(binary_vars), all_of(numeric_vars), "time", "status")

# Convert status variable from 1 and 2 to 0 and 1 (in the lung dataset, 1 indicates death and 2 indicates survival)
data1$status <- ifelse(data1$status == 1, 1, 0)

# Remove missing values
data1 <- na.omit(data1)

# Print the structure of the dataframe
str(data1)

# Check unique values of binary variables
sapply(data1[binary_vars], unique)

# Convert sex variable from "1" to 1 and "2" to 0 (in the lung dataset, 1 indicates male and 2 indicates female)
data1$sex <- ifelse(data1$sex == 1, 1, 0)

# Confirm the conversion results
str(data1[binary_vars])

# Create survival object
surv_obj <- Surv(data1$time, data1$status)
print(surv_obj)

# Ensure reproducibility by setting the seed
set.seed(123)

# Run random forest survival model
rf_model <- ranger(
  formula = Surv(time, status) ~ .,  # Survival model formula
  data = data1,                      # Data
  importance = 'permutation',        # Variable importance calculation method
  num.trees = 100,                   # Number of trees
  seed = 123                         # Random seed
)

# Print model results
print(rf_model)

# Get variable importance
var_importance <- rf_model$variable.importance

# Print the top 20 variables by importance
var_importance_sorted <- sort(var_importance, decreasing = TRUE)
print(var_importance_sorted[1:20])

# Select the top 20 important variables
top_vars <- names(var_importance_sorted[1:20])
print(top_vars)

# Convert variable importance to dataframe
importance_df <- data.frame(
  Variable = names(var_importance),
  Importance = var_importance
)

# Sort the importance and select the top 20 variables
importance_df <- importance_df %>%
  arrange(desc(Importance)) %>%
  head(20)  # Select the top 20 variables

# Plot the top 20 variables by importance
ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip the coordinates so the variable names are on the y-axis for better readability
  labs(
    title = "Variable Importance from Random Forest Model",
    x = "Variables",
    y = "Importance"
  ) +
  theme_minimal()








# method 2--using randomForestSRC and satpred

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
