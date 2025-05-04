
######################
#     FINANCE 4      #
#   Homework 1    #
#  Yaru(Janet) YANG   #

# Install only if not already installed
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("dbplyr")) install.packages("dbplyr")
if (!require("broom")) install.packages("broom")

##### LIBRARIES #####
set.seed(12345) #setup
library(tidyverse)
library(broom)
library(ggplot2)

#Omitted variables bias####
#simulation data
data = function(n = 10000, beta_zx = 1, beta_yx = 1, beta_yz = 2) {
  x = rnorm(n)
  eta = rnorm(n)
  z = beta_zx * x + eta
  epsilon = rnorm(n)
  y = 1 + beta_yx * x + beta_yz * z + epsilon
  
  true_mod = lm(y ~ x + z)
  omit_mod = lm(y ~ x)
  
  tibble(
    beta_yx_true = coef(true_mod)["x"],
    beta_yx_omitted = coef(omit_mod)["x"]
  )
}
#repeated simulations fpr 1000 times
results = replicate(1000, data(), simplify = FALSE) %>% bind_rows()
#collect main statistics
summary_stats = results %>%
  summarise(
    mean_true = mean(beta_yx_true),
    mean_omitted = mean(beta_yx_omitted),
    bias = mean(beta_yx_omitted - beta_yx_true),
    sd_omitted = sd(beta_yx_omitted)
  )
#table
knitr::kable(summary_stats, digits = 3)
knitr::kable(summary_stats, digits = 3, format = "latex")


#show the convergence
ggplot(results, aes(x = beta_yx_omitted)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 30) +
  geom_vline(aes(xintercept = mean(beta_yx_omitted)), color = "red", linetype = "dashed") +
  labs(title = "Distribution of Coefficient on x (Omitted z)", x = expression(hat(beta)[yx]^omitted))
test_params = expand_grid(
  beta_zx = c(0.7, 0.8, 0.9),
  beta_yz = c(0.3, 0.4, 0.5)
)

#sensitivity test
sensitivity = test_params %>%
  mutate(result = map2(beta_zx, beta_yz, ~ {
    replicate(500, data(beta_zx = .x, beta_yz = .y), simplify = FALSE) %>% bind_rows() %>%
      summarise(mean_bias = mean(beta_yx_omitted - beta_yx_true))
  })) %>%
  unnest(result)
#table
knitr::kable(sensitivity, digits = 3)
knitr::kable(sensitivity, digits = 3, format = "latex")
#plot 
ggplot(sensitivity, aes(x = factor(beta_zx), y = factor(beta_yz), fill = mean_bias)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Omitted Variable Bias Across Parameter Settings",
    x = expression(beta[zx]),
    y = expression(beta[yz]),
    fill = "Mean Bias"
  ) +
  theme_minimal(base_size = 14)

#collider bias####
#simulation data
data2 = function(n = 10000, beta_yx = 1, beta_cx = 1, beta_cy = 1) {
  x = rnorm(n)
  epsilon = rnorm(n)
  y = 1 + beta_yx * x + epsilon
  
  nu = rnorm(n)
  c = 1 + beta_cx * x + beta_cy * y + nu
  
  true_mod = lm(y ~ x)
  col_mod = lm(y ~ x + c)
  
  tibble(
    beta_yx_true = coef(true_mod)["x"],
    beta_yx_collider = coef(col_mod)["x"]
  )
}

# Repeat the simulation 
results2 = replicate(1000, data2(), simplify = FALSE) %>% bind_rows()
#collect main statistics
summary_stats_2 = results2 %>%
  summarise(
    mean_true = mean(beta_yx_true),
    mean_collider = mean(beta_yx_collider),
    bias = mean(beta_yx_collider - beta_yx_true),
    sd_collider = sd(beta_yx_collider)
  )

# table
knitr::kable(summary_stats_2, digits = 3, caption = "Collider Bias: Mis-specified vs. Correct Model")
knitr::kable(summary_stats_2, digits = 5, format = "latex")
# Plot
results_long = results2 %>%
  pivot_longer(everything(), names_to = "Model", values_to = "Estimate")


#survivor bias####
#import data from WRDS

ff = read_csv("https://raw.githubusercontent.com/Janetyaruyang/endogeneity-simulation/refs/heads/main/FamaFrench.csv") %>%
  mutate(
    rm = mktrf + rf
  ) %>%
  tail(120)


#simulation function
data3 = function(mktrf, rm) {
  n_funds = 1000
  n_months = length(mktrf)
  
  # set beta and eps
  betas = rnorm(n_funds, mean = 1, sd = 0.25)
  eps = matrix(rnorm(n_funds * n_months, mean = 0, sd = 0.01), nrow = n_months, ncol = n_funds)
  
  # get R_it
  factor_matrix = outer(mktrf, betas, "*")  # matrix: 120 × 1000
  R_it = factor_matrix + eps  # set alpha as zero
  
  # survivor setting
  alive = matrix(TRUE, nrow = n_months, ncol = n_funds)
  dead_status = rep(FALSE, n_funds)
  
  # decide the label
  for (year in 0:9) {
    #ten years
    start = year * 12 + 1
    end = start + 11
    fund_yearly = colSums(R_it[start:end, ]) #fund return
    market_yearly = sum(rm[start:end]) #market return
    
    dead_this_year = fund_yearly < (market_yearly - 0.10)
    dead_status = dead_status | dead_this_year
    
    if (end < n_months) {
      alive[(end + 1):n_months, dead_status] = FALSE
    }
  }
  
  # regression to get alpha
  R_it[!alive] = NA
  alphas = map_dbl(1:n_funds, function(i) {
    if (sum(!is.na(R_it[, i])) < 12) return(NA)
    fit = lm(R_it[, i] ~ mktrf)
    coef(fit)[1]
  })
  
  # output
  tibble(
    mean_alpha_survivors = mean(alphas[!dead_status], na.rm = TRUE),
    mean_alpha_dead = mean(alphas[dead_status], na.rm = TRUE),
    mean_alpha_all = mean(alphas, na.rm = TRUE)
  )
}

#simulation
set.seed(123) 

mktrf = ff$mktrf
rm = ff$rm

results3 = replicate(1000, data3(mktrf, rm), simplify = FALSE) %>%
  bind_rows()

summary_stats3 = results3 %>%
  summarise(
    avg_survivors = mean(mean_alpha_survivors, na.rm = TRUE),
    avg_dead = mean(mean_alpha_dead, na.rm = TRUE),
    avg_all = mean(mean_alpha_all, na.rm = TRUE),
    sd_all = sd(mean_alpha_all, na.rm = TRUE)
  )

knitr::kable(summary_stats3, digits = 6,caption = "Average Alpha Across 1000 Simulations")
knitr::kable(summary_stats3, digits = 6, format = "latex")
