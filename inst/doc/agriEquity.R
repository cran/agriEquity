## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)

## -----------------------------------------------------------------------------
library(agriEquity)

## -----------------------------------------------------------------------------
bi <- benefit_incidence(farm_subsidy, "fertiliser_subsidy",
                        group = "size_class", weights = "weight",
                        area = "land_ha", income = "farm_income")
bi[, c("group", "population_share", "benefit_share", "parity_ratio",
       "benefit_per_hectare")]

## -----------------------------------------------------------------------------
with(farm_subsidy, progressivity_index(fertiliser_subsidy, farm_income,
                                       weights = weight))
with(farm_subsidy, progressivity_index(income_transfer, farm_income,
                                       weights = weight))

## -----------------------------------------------------------------------------
se <- subsidy_equity(farm_subsidy, "fertiliser_subsidy",
                     dimensions = c("size_class", "gender", "region"),
                     weights = "weight",
                     target = list(size_class = c("marginal", "small"),
                                   gender = "female",
                                   region = c("Rainfed-East", "Hill-Tribal")))
se$summary

## -----------------------------------------------------------------------------
de <- with(farm_subsidy, distributional_effects(
  farm_income,
  list(fertiliser = fertiliser_subsidy, transfer = income_transfer,
       power = power_subsidy),
  weights = weight, poverty_line = 60000))
de$summary[, c("scenario", "gini", "fgt0", "fgt1", "lifted_out_share",
               "cost_per_household_lifted")]

## -----------------------------------------------------------------------------
fs <- farm_subsidy
prio <- with(fs, 1 + (size_class %in% c("marginal", "small")) +
               0.5 * (gender == "female") +
               0.5 * (region %in% c("Rainfed-East", "Hill-Tribal")))
br <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
                          weights = fs$weight, area = fs$land_ha,
                          priority = prio,
                          target = fs$size_class %in% c("marginal", "small"),
                          return_rate = fs$marginal_return,
                          env_cost = fs$env_cost_per_rupee,
                          inequality_aversion = 1, env_weight = 1,
                          poverty_line = 60000)
m <- br$metrics
m[, c("rule", "productivity_gain", "environmental_cost", "gini_post",
      "target_share", "poverty_headcount")]

## -----------------------------------------------------------------------------
tf <- tradeoff_frontier(fs$farm_income, fs$fertiliser_subsidy,
                        weights = fs$weight,
                        target = fs$size_class %in% c("marginal", "small"),
                        return_rate = fs$marginal_return,
                        env_cost = fs$env_cost_per_rupee,
                        aversion = c(0, 0.5, 1, 2), env_weights = c(0, 1, 3))
tf
plot(tf)

