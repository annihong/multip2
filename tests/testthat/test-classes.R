
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
