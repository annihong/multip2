library(testthat)
#source("R/helpers.R")

# Test dimension_check function
test_that("dimension_check returns TRUE when all matrices have the same dimensions", {
  x <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
  n = 2
  expect_true(dimension_check(x, n))
})

test_that("dimension_check returns FALSE when matrices have different dimensions", {
  x <- list(matrix(1:4, nrow = 2), matrix(5:9, nrow = 3))
   n = 3
  expect_false(dimension_check(x, n))
})

# Test binary_check function
test_that("binary_check returns TRUE when all matrices are binary", {
  x <- list(matrix(c(0, 1, 1, 0), nrow = 2), matrix(c(1, 0, 0, 1), nrow = 2))
  expect_true(binary_check(x))
})

test_that("binary_check returns FALSE when matrices are not binary", {
  x <- list(matrix(c(0, 1, 1, 2), nrow = 2), matrix(c(1, 0, 0, 1), nrow = 2))
  expect_false(binary_check(x))
})

# Test na_check function
test_that("na_check returns TRUE when there are no missing values", {
  x <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
  expect_true(na_check(x))
})

test_that("na_check returns FALSE when there are missing values", {
  x <- list(matrix(c(1, NA, 3, 4), nrow = 2), matrix(5:8, nrow = 2))
  expect_false(na_check(x))
})

# Test is_valid function
test_that("is_valid returns TRUE for dep_net type with valid input", {
  x <- list(matrix(c(0, 1, 1, 0), nrow = 2), matrix(c(1, 0, 0, 1), nrow = 2))
  n = 2
  expect_true(is_valid(x, "dep_net", n))
})

test_that("is_valid returns FALSE for dep_net type with invalid input", {
  x <- list(matrix(c(0, 1, 1, 2), nrow = 2), matrix(c(1, 0, 0, 1), nrow = 2))
  n = 2
  expect_false(is_valid(x, "dep_net", 2))
})

test_that("is_valid returns TRUE for dyad_covar type with valid input", {
  x <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
  n = 2
  expect_true(is_valid(x, "dyad_covar",n = 2))
})

test_that("is_valid returns FALSE for dyad_covar type with invalid input", {
  x <- list(matrix(c(1, NA, 3, 4), nrow = 2), matrix(5:8, nrow = 2))
  n = 2
  expect_false(is_valid(x, "dyad_covar", n = 2))
})

test_that("is_valid returns TRUE for actor_covar type with valid input", {
  x <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  n = 3
  expect_true(is_valid(x, "actor_covar", 3))
})

test_that("is_valid returns FALSE for actor_covar type with invalid input", {
  x <- data.frame(a = c(1, NA, 3), b = c(4, 5, 6))
  n = 3
  expect_false(is_valid(x, "actor_covar", 3))
})

#tests for matrix_to_dyads_1d
test_that("matrix_to_dyads_1d handles NA values correctly", {
  M <- matrix(c(1, NA, 1, 1), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), c(NA))
})

test_that("matrix_to_dyads_1d returns correct outcomes", {
  M <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), 1)
  
  M <- matrix(c(1, 0, 1, 0), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), 2)
  
  M <- matrix(c(0, 1, 1, 0), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), 4)
  
  M <- matrix(c(0, 0, 0, 0), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), 1)

  M <- matrix(c(0, 1, 0, 0), nrow = 2)
  expect_equal(matrix_to_dyads_1d(M), 3)
})

test_that("matrix_to_dyads_1d returns a vector of the correct length", {
  n = 10
  N = n*(n - 1) / 2
  M <- matrix(1, nrow = n, ncol = n)
  expect_equal(length(matrix_to_dyads_1d(M)), N)
  expect_equal(matrix_to_dyads_1d(M), rep(4, N))
})

test_that("matrix_to_dyads_nd handles uniplex case correctly", {
  M <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equal(matrix_to_dyads_nd(list(M)), 1)
})

test_that("matrix_to_dyads_nd handles multiplex case correctly", {
  M1 <- matrix(c(1, 0, 0, 1), nrow = 2)
  M2 <- matrix(c(0, 1, 1, 0), nrow = 2)
  expect_equal(matrix_to_dyads_nd(list(M1, M2)), 13)
  M1 <- matrix(c(1, 1, 0, 1), nrow = 2)
  M2 <- matrix(c(0, 1, 1, 0), nrow = 2)
  expect_equal(matrix_to_dyads_nd(list(M1, M2)), 15)
})

test_that("matrix_to_dyads_nd returns a vector of the correct length", {
  n = 12
  N = n*(n - 1) / 2
  M1 <- matrix(1, nrow = n, ncol = n)
  M2 <- matrix(1, nrow = n, ncol = n)
  expect_equal(length(matrix_to_dyads_nd(list(M1, M2))), N)
})



test_that("matrix_to_dyads_nd works with larger n", {
  n = 5
  N = n*(n - 1) / 2
  M1 <- matrix(1, nrow = n, ncol = n)
  M2 <-matrix(1, nrow = n, ncol = n)

  expect_equal(matrix_to_dyads_nd(list(M1, M2)), rep(16, N))
})

test_that("matrix_to_dyads_nd handles NA", {
  n = 5
  N = n*(n - 1) / 2
  M1 <- matrix(1, nrow = n, ncol = n)
  M2 <-matrix(1, nrow = n, ncol = n)
  M1[2,1] <- NA
  res <- rep(16, N)
  res[1] <- NA
  expect_equal(matrix_to_dyads_nd(list(M1, M2)), res)
})
