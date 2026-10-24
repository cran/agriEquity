test_that("benefit incidence shares and parity are correct", {
  d <- data.frame(g = c("a", "a", "b", "b"), b = c(10, 0, 30, 60),
                  w = c(1, 1, 1, 1), area = c(1, 1, 2, 4),
                  y = c(100, 100, 200, 300))
  bi <- benefit_incidence(d, "b", group = "g", weights = "w", area = "area",
                          income = "y")
  expect_equal(bi$population_share, c(0.5, 0.5))
  expect_equal(bi$benefit_share, c(0.1, 0.9))
  expect_equal(bi$parity_ratio, c(0.2, 1.8))
  expect_equal(bi$coverage, c(0.5, 1))
  expect_equal(bi$mean_benefit_per_beneficiary, c(10, 45))
  expect_equal(bi$benefit_per_hectare, c(5, 15))
  expect_equal(bi$benefit_income_ratio, c(10 / 200, 90 / 500))
  q <- benefit_incidence(farm_subsidy, "income_transfer",
                         rank_var = "farm_income", weights = "weight")
  expect_equal(nrow(q), 5)
  expect_equal(sum(q$benefit_share), 1)
  expect_equal(sum(q$population_share), 1)
  expect_error(benefit_incidence(d, "b"))
  expect_error(benefit_incidence(d, "zz", group = "g"))
})

test_that("Gini and concentration match known values", {
  # unweighted Gini of 1:4 with the covariance formula
  y <- c(1, 2, 3, 4)
  p <- progressivity_index(benefit = y, income = y)
  expect_equal(p$estimates[["gini_income"]], 0.25)
  expect_equal(p$estimates[["concentration"]], 0.25)
  expect_equal(p$estimates[["kakwani"]], 0)
  # weights equal to replication
  p2 <- progressivity_index(c(1, 2), c(1, 2), weights = c(2, 2))
  p3 <- progressivity_index(c(1, 1, 2, 2), c(1, 1, 2, 2))
  expect_equal(p2$estimates[["gini_income"]], p3$estimates[["gini_income"]])
})

test_that("progressivity signs and decomposition behave", {
  inc <- c(10, 20, 30, 40, 50)
  pro <- progressivity_index(c(5, 4, 3, 2, 1), inc)
  reg <- progressivity_index(c(1, 2, 3, 4, 10), inc)
  expect_lt(pro$estimates[["concentration"]], 0)
  expect_gt(pro$estimates[["kakwani"]], 0)
  expect_lt(reg$estimates[["kakwani"]], 0)
  e <- pro$estimates
  expect_equal(e[["redistributive_effect"]],
               e[["reynolds_smolensky"]] - e[["reranking"]])
  set.seed(2)
  pb <- with(farm_subsidy[1:300, ],
             progressivity_index(fertiliser_subsidy, farm_income,
                                 weights = weight, n_boot = 30))
  expect_true(all(pb$ci[, 1] <= pb$ci[, 2]))
  expect_output(print(pb), "regressive")
  expect_error(progressivity_index(1:3, 1:2))
})
