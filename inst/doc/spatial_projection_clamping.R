## ----load, include=FALSE------------------------------------------------------
library(multiScaleR)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  message = FALSE,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library(multiScaleR)
library(terra)
library(MASS)

## ----simulate-landscape-------------------------------------------------------
set.seed(321)

# Simulate a raster landscape with two binary covariates
rs <- sim_rast(user_seed = 999, dim = 500, resolution = 30)
rs <- terra::subset(rs, c("bin1", "bin2"))

# Apply the TRUE smoothing used to generate the ecological signal
true_sigmas <- c(400, 200)

true_smoothed <- kernel_scale.raster(
  raster_stack = rs,
  sigma = true_sigmas,
  kernel = "gaussian",
  scale_center = FALSE,
  verbose = FALSE
)

## ----standardize-truth--------------------------------------------------------
z_stack <- scale(true_smoothed)

## ----define-truth-------------------------------------------------------------
# Define true model coefficients
alpha <- 0.5
beta1 <- 1.25
beta2 <- -0.75

# Calculate the linear predictor
linpred_true <- alpha + (beta1 * z_stack$bin1) + (beta2 * z_stack$bin2)

# Saturate the linear predictor so extreme gradients remain interpretable
linpred_true <- 4 * tanh(linpred_true / 4)

# Convert to the expected count surface (Poisson mean)
true_mu <- exp(linpred_true)

## ----restrict-sampling--------------------------------------------------------
# Define a restricted environmental envelope for sampling
q1 <- quantile(values(z_stack$bin1), probs = c(0.50, 0.90), na.rm = TRUE)
q2 <- quantile(values(z_stack$bin2), probs = c(0.25, 0.75), na.rm = TRUE)

# Create a mask that isolates this specific domain
sample_mask <- z_stack$bin1
sample_mask[] <- ifelse(
  z_stack$bin1[] >= q1[1] & z_stack$bin1[] <= q1[2] &
    z_stack$bin2[] >= q2[1] & z_stack$bin2[] <= q2[2],
  1, NA
)

# Sample only within this restricted domain
pts <- spatSample(
  sample_mask,
  size = 150,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE
)

# Extract the true mean surface at sampled points and simulate counts
mu_pts <- terra::extract(true_mu, pts)[, 2]
counts <- rpois(length(mu_pts), lambda = mu_pts)

## ----visualize-sampling-domain, fig.width=9, fig.height=4.5-------------------
par(mfrow = c(1, 2))

# Plot the restricted sampling mask and points
plot(sample_mask, main = "Restricted Sampling Domain")
points(pts, pch = 16, cex = 0.45)

# Plot the true mean surface and points
plot(true_mu, main = "True Mean Surface")
points(pts, pch = 16, cex = 0.35)

par(mfrow = c(1, 1))

## ----fit-model----------------------------------------------------------------
# Prepare multiscale data objects based on the sampled points
kernel_inputs <- kernel_prep(
  pts = pts,
  raster_stack = rs,
  max_D = 1500,
  kernel = "gaussian",
  verbose = FALSE
)

# Combine count data with kernel predictors
dat <- data.frame(counts = counts, kernel_inputs$kernel_dat)

# Fit the base GLM
mod <- glm(counts ~ bin1 + bin2, family = poisson(), data = dat)

## ----model-summary------------------------------------------------------------
# Optimize scales of effect
opt_mod <- multiScale_optim(
  fitted_mod = mod,
  kernel_inputs = kernel_inputs
)

summary(opt_mod)

## ----project-unclamped--------------------------------------------------------
# Project covariates without clamping
r_unclamped <- kernel_scale.raster(
  raster_stack = rs,
  multiScaleR = opt_mod,
  scale_center = TRUE,
  clamp = FALSE,       # Clamping disabled
  verbose = FALSE
)

# Predict the expected counts
pred_unclamped <- predict(r_unclamped, opt_mod$opt_mod, type = "response")

# Calculate error (Predicted - True)
unclamped_error <- pred_unclamped - true_mu

## ----compare-ranges-unclamped-------------------------------------------------
# Extract training values and projected raster values
train_vals <- opt_mod$opt_mod$model[, c("bin1", "bin2")]
proj_vals_unclamped <- as.data.frame(values(r_unclamped))

# Construct a table comparing the ranges
range_tab <- data.frame(
  covariate = c("bin1", "bin2"),
  train_min = c(min(train_vals$bin1, na.rm = TRUE), min(train_vals$bin2, na.rm = TRUE)),
  train_max = c(max(train_vals$bin1, na.rm = TRUE), max(train_vals$bin2, na.rm = TRUE)),
  raster_min_unclamped = c(min(proj_vals_unclamped$bin1, na.rm = TRUE),
                           min(proj_vals_unclamped$bin2, na.rm = TRUE)),
  raster_max_unclamped = c(max(proj_vals_unclamped$bin1, na.rm = TRUE),
                           max(proj_vals_unclamped$bin2, na.rm = TRUE))
)

range_tab

## ----visualize-unclamped, fig.width=9, fig.height=4.5-------------------------
par(mfrow = c(1, 3))
plot(true_mu, main = "True Mean Surface")
plot(pred_unclamped, main = "Unclamped Prediction")
plot(unclamped_error, main = "Unclamped Error")
par(mfrow = c(1, 1))

## ----project-clamped----------------------------------------------------------
# Contracted bounds
r_cn20 <- kernel_scale.raster(
  raster_stack = rs,
  multiScaleR = opt_mod,
  scale_center = TRUE,
  clamp = TRUE,
  pct_mx = -0.20,
  verbose = FALSE
)

# Strict clamping
r_c0 <- kernel_scale.raster(
  raster_stack = rs,
  multiScaleR = opt_mod,
  scale_center = TRUE,
  clamp = TRUE,
  pct_mx = 0,
  verbose = FALSE
)

# Expanded bounds
r_cp20 <- kernel_scale.raster(
  raster_stack = rs,
  multiScaleR = opt_mod,
  scale_center = TRUE,
  clamp = TRUE,
  pct_mx = 0.20,
  verbose = FALSE
)

# Generate predictions for each clamped raster
pred_cn20 <- predict(r_cn20, opt_mod$opt_mod, type = "response")
pred_c0   <- predict(r_c0,   opt_mod$opt_mod, type = "response")
pred_cp20 <- predict(r_cp20, opt_mod$opt_mod, type = "response")

# Calculate errors
error_cn20 <- pred_cn20 - true_mu
error_c0   <- pred_c0   - true_mu
error_cp20 <- pred_cp20 - true_mu

## ----compare-ranges-clamped---------------------------------------------------
proj_vals_cn20 <- as.data.frame(values(r_cn20))
proj_vals_c0   <- as.data.frame(values(r_c0))
proj_vals_cp20 <- as.data.frame(values(r_cp20))

range_tab_clamp <- data.frame(
  covariate = rep(c("bin1", "bin2"), each = 4),
  setting = rep(c("unclamped", "pct_mx = -0.20", "pct_mx = 0", "pct_mx = 0.20"), times = 2),
  min_value = c(
    min(proj_vals_unclamped$bin1, na.rm = TRUE),
    min(proj_vals_cn20$bin1, na.rm = TRUE),
    min(proj_vals_c0$bin1, na.rm = TRUE),
    min(proj_vals_cp20$bin1, na.rm = TRUE),
    min(proj_vals_unclamped$bin2, na.rm = TRUE),
    min(proj_vals_cn20$bin2, na.rm = TRUE),
    min(proj_vals_c0$bin2, na.rm = TRUE),
    min(proj_vals_cp20$bin2, na.rm = TRUE)
  ),
  max_value = c(
    max(proj_vals_unclamped$bin1, na.rm = TRUE),
    max(proj_vals_cn20$bin1, na.rm = TRUE),
    max(proj_vals_c0$bin1, na.rm = TRUE),
    max(proj_vals_cp20$bin1, na.rm = TRUE),
    max(proj_vals_unclamped$bin2, na.rm = TRUE),
    max(proj_vals_cn20$bin2, na.rm = TRUE),
    max(proj_vals_c0$bin2, na.rm = TRUE),
    max(proj_vals_cp20$bin2, na.rm = TRUE)
  )
)

range_tab_clamp

## ----visualize-clamped-predictions, fig.width=10, fig.height=8----------------
par(mfrow = c(2, 2))
plot(pred_unclamped, main = "Unclamped")
plot(pred_cn20, main = "Clamped: pct_mx = -0.20")
plot(pred_c0, main = "Clamped: pct_mx = 0")
plot(pred_cp20, main = "Clamped: pct_mx = 0.20")
par(mfrow = c(1, 1))

## ----visualize-clamped-errors, fig.width=10, fig.height=8---------------------
par(mfrow = c(2, 2))
plot(unclamped_error, main = "Error: Unclamped")
plot(error_cn20, main = "Error: pct_mx = -0.20")
plot(error_c0, main = "Error: pct_mx = 0")
plot(error_cp20, main = "Error: pct_mx = 0.20")
par(mfrow = c(1, 1))

## ----rmse-comparison----------------------------------------------------------
# Create masks representing areas inside and outside the sampled domain
inside_mask <- sample_mask
outside_mask <- sample_mask
inside_mask[] <- ifelse(!is.na(sample_mask[]), 1, NA)
outside_mask[] <- ifelse(is.na(sample_mask[]), 1, NA)

# Helper function to compute RMSE
rmse <- function(pred, truth, mask) {
  p <- values(pred, mat = FALSE)
  t <- values(truth, mat = FALSE)
  m <- values(mask, mat = FALSE)
  
  idx <- !is.na(m) & is.finite(p) & is.finite(t)
  sqrt(mean((p[idx] - t[idx])^2, na.rm = TRUE))
}

# Compile RMSE scores into a table
rmse_tab <- data.frame(
  model = c("unclamped", "pct_mx = -0.20", "pct_mx = 0", "pct_mx = 0.20"),
  RMSE_inside = c(
    rmse(pred_unclamped, true_mu, inside_mask),
    rmse(pred_cn20, true_mu, inside_mask),
    rmse(pred_c0, true_mu, inside_mask),
    rmse(pred_cp20, true_mu, inside_mask)
  ),
  RMSE_outside = c(
    rmse(pred_unclamped, true_mu, outside_mask),
    rmse(pred_cn20, true_mu, outside_mask),
    rmse(pred_c0, true_mu, outside_mask),
    rmse(pred_cp20, true_mu, outside_mask)
  )
)

rmse_tab

