setwd("~/Documents/Research/UM_NSF/MatTraits/Maternal-physiology-main")


### EP -------
EP_maternal <- read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx") %>%
  as.data.frame()

EP_maternal <- EP_maternal[-which(EP_maternal$Strain=="LL"),]

columns_to_resid_Lean <- colnames(EP_maternal)[c(12:14,16,18:21)]

EP_maternal_resid <- EP_maternal %>%
  group_by(Strain) %>%
  mutate(across(
    all_of(columns_to_resid_Lean),
    ~ {
      complete_cases <- !is.na(.) & !is.na(Start_AvgLean)
      
      fit <- lm(.[complete_cases] ~ Start_AvgLean[complete_cases])
      
      # Create NA-filled vector of same length as group
      res <- rep(NA_real_, n())
      
      # Assign residuals where data is complete
      res[complete_cases] <- resid(fit)
      res
    },
    .names = "resid_{.col}"
  ),
  resid_Lung_H20 = {
    complete_cases <- !is.na(Lung_H20) & !is.na(Lung_drymass)
    fit <- lm(Lung_H20[complete_cases] ~ Lung_drymass[complete_cases])
    res <- rep(NA_real_, n())
    res[complete_cases] <- resid(fit)
    res
  }
  ) %>%
  ungroup()

columns_to_scale <- colnames(EP_maternal_resid)[12:24]

EP_maternal_scaled <- EP_maternal_resid %>%
                        group_by(Strain) %>%
                        mutate(across(all_of(columns_to_scale), ~ as.numeric(scale(.)), 
                                      .names = "scaled_{.col}")) %>%
                        ungroup()

write.csv(EP_maternal_scaled, "MetaData_EP_Maternal_resid_scaled.csv")




### LP -------
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx") %>%
  as.data.frame()

columns_to_resid_Lean <- colnames(LP_maternal)[c(18:21,24:27)]

LP_maternal_resid <- LP_maternal %>%
  group_by(Strain) %>%
  mutate(across(
    all_of(columns_to_resid_Lean),
    ~ {
      complete_cases <- !is.na(.) & !is.na(Start_AvgLean)
      
      fit <- lm(.[complete_cases] ~ Start_AvgLean[complete_cases])
      
      # Create NA-filled vector of same length as group
      res <- rep(NA_real_, n())
      
      # Assign residuals where data is complete
      res[complete_cases] <- resid(fit)
      res
    },
    .names = "resid_{.col}"
  ),
  resid_Lung_H20 = {
    complete_cases <- !is.na(Lung_H20) & !is.na(Lung_drymass)
    fit <- lm(Lung_H20[complete_cases] ~ Lung_drymass[complete_cases])
    res <- rep(NA_real_, n())
    res[complete_cases] <- resid(fit)
    res
  }
  ) %>%
  ungroup()

columns_to_scale <- colnames(LP_maternal_resid)[18:30]

LP_maternal_scaled <- LP_maternal_resid %>%
  group_by(Strain) %>%
  mutate(across(all_of(columns_to_scale), ~ as.numeric(scale(.)), 
                .names = "scaled_{.col}")) %>%
  ungroup()

write.csv(LP_maternal_scaled, "MetaData_LP_Maternal_resid_scaled.csv")



##EP trait cor
d <- read.csv("Raw_Data/MetaData_EP_Maternal_resid_scaled.csv")

library(dplyr)
library(tidyr)
library(purrr)

# Example: specify the columns to correlate
columns_to_correlate <- c("resid_FoodCon", "resid_F.Mass_Gain", "resid_Lung_drymass", "resid_RV", "resid_LV_Sep", "resid_Lung_H20",  "scaled_Start_AvgLean", "scaled_Hb", "scaled_HctAv", "scaled_Gluc")

# Group by 'group' and calculate correlation matrix for each group
cor_with_p <- function(data) {
  vars <- colnames(data)
  n <- length(vars)
  cor_mat <- matrix(NA, n, n)
  p_mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "pairwise.complete.obs")
      cor_mat[i, j] <- test$estimate
      p_mat[i, j] <- test$p.value
    }
  }
  colnames(cor_mat) <- rownames(cor_mat) <- vars
  colnames(p_mat) <- rownames(p_mat) <- vars
  list(cor = cor_mat, p = p_mat)
}

correlation_by_group <- d %>%
  group_by(Strain) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$Strain))) %>%
  map(~ {
    data <- select(.x, all_of(columns_to_correlate))
    cor_with_p(data)
  })

dir.create("correlation_output", showWarnings = FALSE)

# Loop over groups and write to CSV
walk(names(correlation_by_group), function(group_name) {
  cor_mat <- correlation_by_group[[group_name]]$cor
  p_mat <- correlation_by_group[[group_name]]$p
  
  write.csv(cor_mat, file = paste0("correlation_output/", group_name, "_correlation.csv"), row.names = TRUE)
  write.csv(p_mat, file = paste0("correlation_output/", group_name, "_pvalues.csv"), row.names = TRUE)
})





##LP trait cor --
d <- read.csv("Raw_Data/MetaData_LP_Maternal_resid_scaled.csv")

library(dplyr)
library(tidyr)
library(purrr)

# Example: specify the columns to correlate
columns_to_correlate <- c("resid_FoodCon", "resid_F.Mass_Gain", "resid_Lung_drymass", "resid_RV", "resid_LV_Sep", "resid_Lung_H20",  "scaled_Start_AvgLean", "scaled_Hb", "scaled_HctAv", "scaled_Gluc")

# Group by 'group' and calculate correlation matrix for each group
cor_with_p <- function(data) {
  vars <- colnames(data)
  n <- length(vars)
  cor_mat <- matrix(NA, n, n)
  p_mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "pairwise.complete.obs")
      cor_mat[i, j] <- test$estimate
      p_mat[i, j] <- test$p.value
    }
  }
  colnames(cor_mat) <- rownames(cor_mat) <- vars
  colnames(p_mat) <- rownames(p_mat) <- vars
  list(cor = cor_mat, p = p_mat)
}

correlation_by_group <- d %>%
  group_by(Strain) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$Strain))) %>%
  map(~ {
    data <- select(.x, all_of(columns_to_correlate))
    cor_with_p(data)
  })

dir.create("LP_correlation_output", showWarnings = FALSE)

# Loop over groups and write to CSV
walk(names(correlation_by_group), function(group_name) {
  cor_mat <- correlation_by_group[[group_name]]$cor
  p_mat <- correlation_by_group[[group_name]]$p
  
  write.csv(cor_mat, file = paste0("LP_correlation_output/", group_name, "_correlation.csv"), row.names = TRUE)
  write.csv(p_mat, file = paste0("LP_correlation_output/", group_name, "_pvalues.csv"), row.names = TRUE)
})


