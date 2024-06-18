test_that("Uniplex empty model works as expected", {
  # Create test data
  dep_net <- list(matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 10, ncol = 10))
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function
  model <- Mp2Model(dep_net, dyad_covar, actor_covar)
  names(dep_net) <- "network1" 

  # Check the attributes of the model
  expect_equal(model$t, 1)
  expect_equal(model$H, 0)
  expect_equal(model$n, 10)
  expect_equal(model$dep_lab, "network1")
  expect_equal(model$pair_lab, NULL)
  expect_equal(model$dyad_covar_lab, names(dyad_covar))
  expect_equal(model$actor_covar_lab, names(actor_covar))
  expect_equal(model$data$dep_net, dep_net)
  expect_equal(model$data$dyad_covar, dyad_covar)
  expect_equal(model$data$actor_covar, actor_covar)
  expect_equal(model$model$covar, list())
  expect_equal(model$model$prior, default_prior_empty_mdl(names(dep_net), NULL))
  expect_equal(model$model$full_prior, NULL)
  expect_equal(model$fit_res$stan_data, NULL)
  expect_equal(model$fit_res$stan_fit, NULL)
  expect_equal(model$fit_res$summary, NULL)
  expect_equal(model$fit_res$par_labels, NULL)
})

test_that("Multiplex empty model works as expected", {
  # Create test data
  dep_net <- list("net1" = matrix(0, nrow = 10, ncol = 10), "net2"=matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 10, ncol = 10))
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function
  model <- Mp2Model(dep_net, dyad_covar, actor_covar)

  # Check the attributes of the model
  expect_equal(model$t, 2)
  expect_equal(model$H, 1)
  expect_equal(model$n, 10)
  expect_equal(model$dep_lab, names(dep_net))
  expect_equal(model$pair_lab, get_pair_names(names(dep_net)))
  expect_equal(model$dyad_covar_lab, names(dyad_covar))
  expect_equal(model$actor_covar_lab, names(actor_covar))
  expect_equal(model$data$dep_net, dep_net)
  expect_equal(model$data$dyad_covar, dyad_covar)
  expect_equal(model$data$actor_covar, actor_covar)
  expect_equal(model$model$covar, list())
  expect_equal(model$model$prior, default_prior_empty_mdl(names(dep_net), get_pair_names(names(dep_net))))
  expect_equal(model$model$full_prior, NULL)
  expect_equal(model$fit_res$stan_data, NULL)
  expect_equal(model$fit_res$stan_fit, NULL)
  expect_equal(model$fit_res$summary, NULL)
  expect_equal(model$fit_res$par_labels, NULL)
})


test_that("Mp2Model throws an error with invalid input", {
  # Create test data with invalid input
  dep_net <- list(matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 5, ncol = 5))  # Invalid dimensions
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function and expect an error
  expect_error(Mp2Model(dep_net, dyad_covar, actor_covar), "is_valid")
})
# # Test for baseline parameter with mean and sd
# test_that("update_prior sets baseline parameter with mean and sd", {
#   model_obj <- list(
#     model = list(
#       prior = list(
#         mu_mean_prior = c(1, 2, 3),
#         mu_sd_prior = c(0.1, 0.2, 0.3)
#       )
#     ),
#     dep_lab = c("density", "reciprocity")
#   )
  
#   updated_model_obj <- update_prior(model_obj, "density", "baseline", mean = 10, sd = 0.5)
  
#   expect_equal(updated_model_obj$model$prior$mu_mean_prior, c(10, 2, 3))
#   expect_equal(updated_model_obj$model$prior$mu_sd_prior, c(0.5, 0.2, 0.3))
# })

# # Test for random parameter with eta, alpha, and beta
# test_that("update_prior sets random parameter with eta, alpha, and beta", {
#   model_obj <- list(
#     model = list(
#       prior = list(
#         LJK_eta_prior = 0.1,
#         scale_alpha_prior = 0.2,
#         scale_beta_prior = 0.3
#       )
#     )
#   )
  
#   updated_model_obj <- update_prior(model_obj, "density", "random", eta = 0.5, alpha = 0.8, beta = 0.9)
  
#   expect_equal(updated_model_obj$model$prior$LJK_eta_prior, 0.5)
#   expect_equal(updated_model_obj$model$prior$scale_alpha_prior, 0.8)
#   expect_equal(updated_model_obj$model$prior$scale_beta_prior, 0.9)
# })

# # Test for covariate parameter with mean and sd
# test_that("update_prior sets covariate parameter with mean and sd", {
#   model_obj <- list(
#     model = list(
#       prior = list(),
#       covar = list(
#         density_covar_lab_layer_lab = list(
#           mean_prior = c(1, 2, 3),
#           sd_prior = c(0.1, 0.2, 0.3)
#         )
#       )
#     )
#   )
  
#   updated_model_obj <- update_prior(model_obj, "density", "covariate", covar_lab = "covar_lab", mean = 10, sd = 0.5)
  
#   expect_equal(updated_model_obj$model$covar$density_covar_lab_layer_lab$mean_prior, c(10, 2, 3))
#   expect_equal(updated_model_obj$model$covar$density_covar_lab_layer_lab$sd_prior, c(0.5, 0.2, 0.3))
# })

# # Test for invalid type
# test_that("update_prior throws an error for invalid type", {
#   model_obj <- list(
#     model = list(
#       prior = list()
#     )
#   )
  
#   expect_error(update_prior(model_obj, "density", "invalid_type"))
# })

# # Test for missing covariate label
# test_that("update_prior throws an error for missing covariate label", {
#   model_obj <- list(
#     model = list(
#       prior = list()
#     )
#   )
  
#   expect_error(update_prior(model_obj, "density", "covariate"))
# })
