LASSO/岭回归（glmnet）
library(shellpipes)
library(ggplot2)  
library(glmnet)
library(survival)
library(readxl)
library(dplyr)

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

# =
x <- model.matrix(~ . - survival_time - censor, data = train_df)[, -1]
y <- Surv(train_df$survival_time, train_df$censor)

# Use cross-validation to select the best lambda
cv_glmnet_fit <- cv.glmnet(x, y, family = "cox")

# Best lambda
best_lambda <- cv_glmnet_fit$lambda.min
cat("Best lambda: ", best_lambda, "\n")

# Fit glmnet model with best lambda
glmnet_fit <- glmnet(x, y, family = "cox", lambda = best_lambda)
glmnet_coef <- coef(glmnet_fit)
glmnet_coef <- data.frame(coef = rownames(glmnet_coef), value = as.vector(glmnet_coef), model = "glmnet")

# Fit coxph model
coxph_fit <- coxph(y ~ x, method = "breslow")
coxph_coef <- coef(coxph_fit)
coxph_coef <- data.frame(coef = names(coxph_coef), value = as.vector(coxph_coef), model = "coxph")

# Combine coefficient data
coef_df <- do.call("rbind", list(glmnet_coef, coxph_coef))
coef_df$coef <- gsub("V", "x", coef_df$coef)

# Limit the number of variables displayed
top_n <- 20
top_coef_df <- coef_df %>%
  group_by(model) %>%
  slice_max(order_by = abs(value), n = top_n) %>%  
  ungroup()

# Plot coefficients
p1 <- ggplot(top_coef_df, aes(x = reorder(coef, -value), y = value, col = model)) +
  geom_point(aes(shape = model), alpha = 0.2) +
  scale_colour_manual(breaks = c("coxph", "glmnet"), 
                      values = c("coxph" = "black", "glmnet" = "blue")) +
  scale_shape(guide = FALSE) +
  labs(colour = "Model") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6, angle = 0, hjust = 1))

print(p1)
top_20_variable_names <- unique(top_coef_df$coef)
print("Top 20 variables by absolute value of coefficients:")
print(top_20_variable_names)

# Create a table for the top 20 variables
top_20_variable_names <- unique(top_coef_df$coef)
top_20_table <- top_coef_df %>% filter(coef %in% top_20_variable_names)

print("Top 20 variables by absolute value of coefficients:")
print(top_20_variable_names)

# Display the top 20 table
top_20_table <- top_20_table %>%
  arrange(desc(abs(value)))

print(top_20_table)