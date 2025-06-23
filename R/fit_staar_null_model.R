library(ordinal)
library(stats)

# #########################################################################
# The main user-facing function
# #########################################################################
#' #' Fit a Null Model for STAAR, with Support for Ordinal Outcomes
#'
#' This function provides a unified interface to fit null models for the STAAR
#' package. It handles standard GLMs for binary/continuous traits and correctly
#' fits and converts ordinal regression models (proportional odds models) for
#' use with STAAR.
#'
#' @param fixed A formula specifying the null model.
#' @param data A data frame containing all variables in the formula.
#' @param family The model family. Use standard `glm` families like `binomial()` or
#'   `gaussian()`, or use the special string `"ordinal"` for proportional odds models.
#' @param ... Additional arguments passed to `ordinal::clm()` or `stats::glm()`.
#'
#' @return A STAAR-compatible null model object.
#'
#' @importFrom stats predict formula model.matrix gaussian binomial glm
#' @importFrom ordinal clm
#'
#' @export
fit_staar_null_model <- function(fixed, data, family = binomial(link = "logit"), ...){
  is_ordinal <- is.character(family) && family == "ordinal"

  if (is_ordinal) {
    clm_fit <- ordinal::clm(formula = fixed, data = data, ...)
    obj_nullmodel <- create_staar_null_from_clm(clm_obj = clm_fit, full_data = data)

  } else {
    obj_nullmodel <- stats::glm(formula = fixed, data = data, family = family, ...)
    obj_nullmodel$relatedness <- FALSE
  }

  return(obj_nullmodel)
}

# #########################################################################
# Internal converter function
# #########################################################################
create_staar_null_from_clm <- function(clm_obj, full_data) {
  # --- 1. Calculate the "Working" Components for the Score Test ---
  pred_probs <- predict(clm_obj, newdata = full_data, type = "prob")$fit
  y_factor <- full_data[[as.character(formula(clm_obj)[[2]])]]
  indicator_matrix <- model.matrix(~ y_factor - 1)
  working_residual <- rowSums(indicator_matrix - pred_probs)

  # Calculate the linear predictor (eta = X*beta)
  X_matrix <- model.matrix(object = formula(clm_obj), data = clm_obj$model)
  X_no_intercept <- X_matrix[, -which(colnames(X_matrix) == "(Intercept)")]
  latent_predictor <- X_no_intercept %*% clm_obj$beta

  # z = eta + working_residual.
  z_working_variable <- latent_predictor + working_residual

  # Calculate the weights. Var(z) is approximately 1/w.
  categories <- 1:ncol(pred_probs)
  expected_y <- pred_probs %*% categories
  expected_y_squared <- pred_probs %*% (categories^2)
  var_y <- expected_y_squared - (expected_y^2)
  model_weights <- 1 / (var_y + 1e-8)

  dispersion_est <- 1.0

  # --- 2. Assemble the Final "Pseudo glm" Object ---
  staar_null_obj <- list()
  staar_null_obj$y <- as.vector(z_working_variable)
  staar_null_obj$fitted.values <- as.vector(latent_predictor)

  staar_null_obj$family <- gaussian(link = "identity")
  staar_null_obj$weights <- as.vector(model_weights) # Use the variance-based weights
  staar_null_obj$relatedness <- FALSE
  staar_null_obj$terms <- clm_obj$terms
  staar_null_obj$model <- clm_obj$model
  staar_null_obj$summary <- list(dispersion = dispersion_est)
  class(staar_null_obj) <- c("glm", "lm")

  # --- 3. Add Ancillary Information ---
  staar_null_obj$residuals <- as.vector(working_residual)
  staar_null_obj$coefficients <- clm_obj$coefficients
  staar_null_obj$formula <- formula(clm_obj)

  return(staar_null_obj)
}




