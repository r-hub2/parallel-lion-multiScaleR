## ----setup, include=FALSE, echo=FALSE-----------------------------------------
library(multiScaleR)
# Define the base URL for the extdata directory
base_url <- "https://raw.githubusercontent.com/wpeterman/multiScaleR/master/inst/extdata/"

f <- c("opt_exp.RData", "opt_expow.RData", "opt_fixed.RData",
       "opt_umf_counts.RData", "opt_umf_occ.RData", "opt1.RData", 
       "opt2.RData", "opt3.RData", "opt4.RData", "opt5.RData", 
       "opt6.RData", "sim_opt.RData", "zinb.RData")

# Download each file
for (file in f) {
  download.file(
    url = paste0(base_url, file),
    destfile = file.path(tempdir(), file),
    mode = "wb"
  )
  }

rdat <- list.files(tempdir(), pattern = "*.RData", full.names = T)
invisible(lapply(rdat, load, .GlobalEnv))

## ----echo=FALSE---------------------------------------------------------------
plot_kernel(prob = 0.9, 
            sigma = 100, 
            kernel = "gaus", 
            scale_dist = F,
            add_label = F)

## ----echo=FALSE---------------------------------------------------------------
plot_kernel(prob = 0.9, 
            sigma = 50, 
            kernel = "exp", 
            scale_dist = F,
            add_label = F)

## ----echo=FALSE---------------------------------------------------------------
plot_kernel(prob = 0.9, 
            sigma = 250, 
            beta = 5,
            kernel = "expow", 
            scale_dist = F,
            add_label = F)

## ----echo=FALSE---------------------------------------------------------------
plot_kernel(prob = 0.9, 
            sigma = 300, 
            kernel = "fixed", 
            scale_dist = F,
            add_label = F)

## -----------------------------------------------------------------------------
library(multiScaleR)

## Read in data
data("landscape_counts")
dat <- landscape_counts

data("surv_pts")
pts <- vect(surv_pts)

land_rast <- terra::rast(system.file("extdata", 
                                     "landscape.tif", 
                                     package = 'multiScaleR'))

## -----------------------------------------------------------------------------
summary(dat)

pts

land_rast

plot(land_rast)

## Plot with points
plot(land_rast$land1) 
plot(pts, add = T, pch = 19)

## -----------------------------------------------------------------------------
kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = land_rast,
                             max_D = 1700,
                             kernel = 'gaussian',
                             verbose = FALSE)
kernel_inputs

## -----------------------------------------------------------------------------
df <- data.frame(dat,
                 kernel_inputs$kernel_dat)
str(df)

## -----------------------------------------------------------------------------
mod0 <- glm(counts ~ site + land1 + land2 + land3,
            family = poisson(),
            data = df)
summary(mod0)

## ----opt1, eval=FALSE---------------------------------------------------------
# opt1 <- multiScale_optim(fitted_mod = mod0,
#                          kernel_inputs = kernel_inputs)

## -----------------------------------------------------------------------------
summary(opt1)

## ----opt2---------------------------------------------------------------------
## New model
mod0_2 <- glm(counts ~ site + land1 + land2,
              family = poisson(),
              data = df)

## ----eval=FALSE---------------------------------------------------------------
# ## Optimize
# opt2 <- multiScale_optim(fitted_mod = mod0_2,
#                          kernel_inputs = kernel_inputs)

## -----------------------------------------------------------------------------
summary(opt2)

## -----------------------------------------------------------------------------
## Kernel function
plot(opt2)

## Kernel function; 99% contribution
plot(opt2, prob = 0.99)

## -----------------------------------------------------------------------------
plot_marginal_effects(opt2)

## ----kernel_raster------------------------------------------------------------
rast_opt <- kernel_scale.raster(raster_stack = land_rast,
                                multiScaleR = opt2)
plot(rast_opt)

## ----other_kernel-------------------------------------------------------------
## Negative Exponential
exp_inputs <- kernel_prep(pts = pts,
                          raster_stack = land_rast,
                          max_D = 1700,
                          kernel = 'exp',
                          verbose = FALSE)

## Exponential Power
expow_inputs <- kernel_prep(pts = pts,
                            raster_stack = land_rast,
                            max_D = 1700,
                            kernel = 'expow',
                            verbose = FALSE)

## Fixed width buffer
fixed_inputs <- kernel_prep(pts = pts,
                            raster_stack = land_rast,
                            max_D = 1700,
                            kernel = 'fixed',
                            verbose = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# opt_exp <- multiScale_optim(fitted_mod = mod0_2,
#                             kernel_inputs = exp_inputs)
# 
# ## Starting values needed
# opt_expow <- multiScale_optim(fitted_mod = mod0_2,
#                               kernel_inputs = expow_inputs,
#                               par = c(500/expow_inputs$unit_conv,
#                                       500/exp_inputs$unit_conv,
#                                       5,10))
# 
# opt_fixed <- multiScale_optim(fitted_mod = mod0_2,
#                               kernel_inputs = fixed_inputs)

## -----------------------------------------------------------------------------
mod_list <- list(opt2,
                 opt_exp,
                 opt_expow,
                 opt_fixed)

## AIC table
aic_tab(mod_list)

## BIC table
bic_tab(mod_list)

## -----------------------------------------------------------------------------
## Landscape only effect
mod0_3 <- glm(counts ~ land1 + land2,
              family = poisson(),
              data = df)

## Landscape 1 only effect
mod0_4 <- glm(counts ~ land1,
              family = poisson(),
              data = df)

## Landscape 2 only effect
mod0_5 <- glm(counts ~ land2,
              family = poisson(),
              data = df)

## Landscape 3 only effect
mod0_6 <- glm(counts ~ land3,
              family = poisson(),
              data = df)

## Site only effect
## No multiScaleR optimization
mod0_7 <- glm(counts ~ site,
              family = poisson(),
              data = df)

## ----eval=FALSE---------------------------------------------------------------
# opt3 <- multiScale_optim(fitted_mod = mod0_3,
#                          kernel_inputs = kernel_inputs)
# opt4 <- multiScale_optim(fitted_mod = mod0_4,
#                          kernel_inputs = kernel_inputs)
# opt5 <- multiScale_optim(fitted_mod = mod0_5,
#                          kernel_inputs = kernel_inputs)
# opt6 <- multiScale_optim(fitted_mod = mod0_6,
#                          kernel_inputs = kernel_inputs)

## -----------------------------------------------------------------------------
mod_list2 <- list(opt1, opt2, opt3, opt4, opt5, opt6, mod0_7)

aic_tab(mod_list2)

bic_tab(mod_list2)

## -----------------------------------------------------------------------------
rs <- sim_rast(user_seed = 999, dim = 250)
plot(rs)

## -----------------------------------------------------------------------------
rs <- terra::subset(rs, c(1,4))
s_dat <- sim_dat_unmarked(alpha = 1,
                          beta = c(0.75,-0.75),
                          kernel = 'gaussian',
                          sigma = c(75, 150),
                          n_points = 75,
                          n_surv = 5,
                          det = 0.5,
                          type = 'count',
                          raster_stack = rs,
                          max_D = 550,
                          user_seed = 555)
plot(s_dat$df$y ~ s_dat$df$bin1)
plot(s_dat$df$y ~ s_dat$df$cont2)

## -----------------------------------------------------------------------------
library(unmarked)
kernel_inputs <- kernel_prep(pts = s_dat$pts,
                             raster_stack = rs,
                             max_D = 550,
                             kernel = 'gaus',
                             verbose = FALSE)

umf <- unmarkedFramePCount(y = s_dat$y,
                           siteCovs = kernel_inputs$kernel_dat)

## Base unmarked model
mod0_umf.p <- unmarked::pcount(~1 ~bin1 + cont2,
                               data = umf, 
                               K = 100)

## ----eval=FALSE---------------------------------------------------------------
# opt_umf.p <- multiScale_optim(fitted_mod = mod0_umf.p,
#                               kernel_inputs = kernel_inputs,
#                               n_cores = 8)

## -----------------------------------------------------------------------------
summary(opt_umf.p)

plogis(opt_umf.p$opt_mod@estimates@estimates$det@estimates[[1]])

## -----------------------------------------------------------------------------
plot(opt_umf.p)
plot_marginal_effects(opt_umf.p)

## -----------------------------------------------------------------------------
## Project model 
rast_scale.center <- kernel_scale.raster(raster_stack = rs,
                                         multiScaleR = opt_umf.p,
                                         scale_center = TRUE,
                                         clamp = TRUE,
                                         pct_mx = 0.00)

plot(rast_scale.center)


ab.mod_pred <- terra::predict(rast_scale.center, 
                             opt_umf.p$opt_mod, 
                             type = 'state')

plot(ab.mod_pred)

## -----------------------------------------------------------------------------
s_dat.occ <- sim_dat_unmarked(alpha = 0.25,
                              beta = c(-0.75,1),
                              kernel = 'gaussian',
                              sigma = c(225, 75),
                              n_points = 250,
                              n_surv = 5,
                              det = 0.5,
                              type = 'occ',
                              raster_stack = rs,
                              max_D = 800,
                              user_seed = 555)

plot(s_dat.occ$df$y ~ s_dat.occ$df$bin1)
plot(s_dat.occ$df$y ~ s_dat.occ$df$cont2)

## -----------------------------------------------------------------------------
kernel_inputs <- kernel_prep(pts = s_dat.occ$pts,
                             raster_stack = rs,
                             max_D = 800,
                             kernel = 'gaus',
                             verbose = FALSE)

## Occupancy frame
umf <- unmarkedFrameOccu(y = s_dat.occ$y,
                         siteCovs = kernel_inputs$kernel_dat)

## Base unmarked model
(mod0_umf.occ <- unmarked::occu(~1 ~bin1 + cont2,
                                data = umf))

## ----eval=FALSE---------------------------------------------------------------
# opt_umf.occ <- multiScale_optim(fitted_mod = mod0_umf.occ,
#                                 kernel_inputs = kernel_inputs,
#                                 n_cores = 8)

## -----------------------------------------------------------------------------
summary(opt_umf.occ)
plogis(opt_umf.occ$opt_mod@estimates@estimates$det@estimates[[1]])

## -----------------------------------------------------------------------------
## Project model 
rast_scale.center <- kernel_scale.raster(raster_stack = rs,
                                         multiScaleR = opt_umf.occ,
                                         scale_center = TRUE,
                                         clamp = T)
plot(rast_scale.center)

occ.mod_pred <- terra::predict(rast_scale.center, 
                             opt_umf.occ$opt_mod, 
                             type = 'state')

plot(occ.mod_pred)

## -----------------------------------------------------------------------------
set.seed(12345)
rs <- sim_rast(user_seed = 999, 
               dim = 500, resolution = 30)
rs <- terra::subset(rs, c(1,3))

## Zero-inflated parameters
zi_alpha <- -0.25
zi_beta <- -1

n_points <- 400
alpha <- 0.5
beta <- 0.75
kernel <- 'gaussian'
sigma <- c(400, # For main effect, bin1
           200) # For zero-inflated, cont1
StDev <- 2  # Negative binomial dispersion

## Generate random points
bnd <- buffer(vect(ext(rs[[1]])), -1000)
pts <- spatSample(x = rs[[1]], 
                  size = n_points, 
                  method = 'random',
                  ext = ext(bnd),
                  replace = F,
                  as.points = T)

kernel_out <- kernel_prep(pts = pts,
                          raster_stack = rs,
                          max_D = 1500,
                          kernel = kernel,
                          sigma = sigma,
                          verbose = FALSE)

## ZINB simulation
zi_prob <- plogis(zi_alpha + zi_beta*kernel_out$kernel_dat$cont1)

# Simulate zero-inflated counts
# 1 = excess zero, 0 = from Poisson
is_zero <- rbinom(length(zi_prob), size = 1, prob = zi_prob)  

obs <- exp(alpha + beta*kernel_out$kernel_dat$bin1)

obs.zinb <- ifelse(is_zero, 0,
                   rnbinom(length(is_zero),
                           mu = obs,
                           size = StDev))

dat <- data.frame(cnt = obs.zinb,
                  kernel_out$kernel_dat)
with(dat, plot(cnt ~ bin1))

## -----------------------------------------------------------------------------
library(pscl)
zinb_mod <- pscl::zeroinfl(cnt ~ bin1 | cont1,
                           dist = 'negbin',
                           data = dat)

kernel_inputs <- kernel_prep(pts = pts,
                             kernel = kernel,
                             max_D = 1500,
                             raster_stack = rs,
                             verbose = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# zinb_opt <- multiScale_optim(zinb_mod,
#                              kernel_inputs,
#                              n_cores = 4)

## -----------------------------------------------------------------------------
summary(zinb_opt)
plot(zinb_opt)
plot_marginal_effects(zinb_opt)

## -----------------------------------------------------------------------------
## With fitted model
kernel_dist(opt2)

## Distance of 95% kernel
kernel_dist(opt2, prob = 0.95)


## -----------------------------------------------------------------------------
kernel_dist(kernel = "gaussian", sigma = 100, prob = 0.9)

kernel_dist(kernel = "exp", sigma = 100, prob = 0.9)

kernel_dist(kernel = "expow", sigma = 100, beta = 5, prob = 0.9)

## -----------------------------------------------------------------------------
plot_kernel(kernel = 'exp',
            sigma = 50)

plot_kernel(kernel = 'expow',
            beta = 5,
            sigma = 50)

plot_kernel(kernel = 'gaussian',
            sigma = 50)

plot_kernel(kernel = 'gaussian',
            sigma = 50)

plot_kernel(kernel = 'gaussian',
            sigma = 50,
            scale_dist = F)

plot_kernel(kernel = 'gaussian',
            sigma = 50,
            add_label = F)

## ----error=FALSE, message=FALSE-----------------------------------------------
r_sim1 <- sim_rast(dim = 100,
                   resolution = 30,
                   user_seed = 555)

r_sim2 <- sim_rast(dim = 100,
                   resolution = 30,
                   autocorr_range1 = 100,
                   autocorr_range2 = 1,
                   sill = 10,
                   user_seed = 555)

plot(r_sim1)
plot(r_sim2)

## -----------------------------------------------------------------------------
s_dat <- sim_dat(alpha = 0.25,
                 beta = c(0.75, -0.75),
                 kernel = 'gaussian',
                 sigma = c(350, 200),
                 type = 'count',
                 n_points = 100,
                 raster_stack = terra::subset(r_sim1, c(1,4)),
                 min_D = 250,
                 user_seed = 999)

## -----------------------------------------------------------------------------
plot(s_dat$df$y ~ s_dat$df$bin1)
plot(s_dat$df$y ~ s_dat$df$cont2)

## -----------------------------------------------------------------------------
sim_mod <- glm(y ~ bin1 + cont2,
               family = 'poisson',
               data = s_dat$df)

kernel_inputs <- kernel_prep(pts = s_dat$pts,
                             raster_stack = r_sim1,
                             max_D = 1500,
                             kernel = 'gaussian',
                             verbose = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# sim_opt <- multiScale_optim(fitted_mod = sim_mod,
#                             kernel_inputs = kernel_inputs)

## -----------------------------------------------------------------------------
summary(sim_opt)

plot(sim_opt)

plot(kernel_scale.raster(raster_stack = r_sim1,
                         multiScaleR = sim_opt))

plot_marginal_effects(sim_opt)

