library(tidyverse)
library(MASS)
library(stats)
library(sf)
library(ggspatial)
library(patchwork)
library(lubridate)
library(rstanarm)
library(bayesplot)
library(loo)

figdir <- "C:/Users/phpuenig/Dropbox/deaths_spatial/figures/"

# ---------------------------------------------------------------------------- #
# READ DATA
# ---------------------------------------------------------------------------- #

# Primary dataset:
dat_tot <- readRDS("dat_tot.rds") %>%
  rename(pop_dens = pop_dens_total,
         pop_dens_scale = pop_dens_total_scale) %>%
  mutate(log_pop_dens = log(pop_dens),
         log_pop_dens_scale = scale(log_pop_dens))

## Shapefiles
regions <- readRDS("shp.rds") %>%
  filter(grepl("E", lad19cd)) # filter to England 

# Hold-out set for theta optimisation
samp <- sample(unique(dat_tot$lad19cd),0.4*n_distinct(dat_tot$lad19cd))
dat_opt <- filter(dat_tot, lad19cd %in% samp)
dat_fit <- filter(dat_tot, !lad19cd %in% samp)


# ---------------------------------------------------------------------------- #
# FUNCTIONS
# ---------------------------------------------------------------------------- #

# Plot map
map_theme <- function () {
  theme_minimal() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank()) 
}

basic_map <- function(sf, fill, rate1e5 = F){
  
  if (rate1e5 == T){
    sf <- mutate(sf, fill = !!sym(fill)*1e5/la_pop)
  }else{sf <- mutate(sf, fill = !!sym(fill))}
  
  p <- ggplot(sf, aes(geometry = geometry, fill = fill)) +
    geom_sf() +
    scale_fill_viridis_c() +
    annotation_scale(location = "br") +
    map_theme() 
  return(p)
}

# Saturation function of density
sat_func <- function(pop_dens, theta = 0.46){
  return((2*theta*pop_dens)/(1 + 2*theta*pop_dens + (1 + 4*theta*pop_dens)^0.5))
}

# Inverse to return density
inv_sat_func <- function(s, theta){
  return(uniroot(function(d) sat_func(d,theta) - s, lower = 0, upper = max(dat_tot$pop_dens)))
}

# Plot model coefficients 
plot_coeffs <- function(mod){
print(plot(mod, pars = names(coefficients(mod))[!names(coefficients(mod)) %in% c("(Intercept)","reciprocal_dispersion")]) + 
        geom_vline(xintercept = 0)
)
}

# ---------------------------------------------------------------------------- #
# MAP DEATHS BY LTLA
# ---------------------------------------------------------------------------- #

# dat_tot %>%
#   full_join(regions) %>%
#   basic_map(fill = "n", rate1e5 = T) +
#   labs(title = "COVID-19-associated mortality by LTLA",
#        subtitle = paste0("Total deaths per 100,000 between week beginning",min(dat$wod)," and ",max(dat$wod)),
#        fill = "") -> map_tot
dat_tot %>%
  full_join(regions) %>%
  basic_map(fill = "SMR") +
  labs(title = "Age-standardised mortality ratio (SMR)",
       fill = "") -> map_SMR
dat_tot %>%
  full_join(regions) %>%
  basic_map(fill = "pop_dens", rate1e5 = F) +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "Population density per km-squared",
       fill = "") -> map_dens

# png(file = "./figures/map_mort.png",height = 400, width = 800)
# map_tot + map_SMR
# dev.off()

png(filename = paste0(figdir, "map_mort_dens.png"), height = 4, width = 8, units = "in", res = 600)
map_SMR + map_dens
dev.off()


# ---------------------------------------------------------------------------- #
# FITTING
# ---------------------------------------------------------------------------- #

reg_prior <- normal(c(0,0.1))

# No effect of population density
fA <- n ~ IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + offset(log(E)) 
A <- stan_glm.nb(fA, data = dat_fit, seed = 12345)
# plot(A, plotfun = "trace")
# summary(A)

prior_summary(A)
# pp_check(A, plotfun = "stat", stat = "mean")
# pp_check(A, plotfun = "dens_overlay")


# Linear effect of population density
fB <- n ~ pop_dens_scale + IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + offset(log(E))
B <- stan_glm.nb(fB, data = dat_fit, seed = 12345)
# plot(B, plotfun = "trace")
# summary(B)

prior_summary(B)
pp_check(B, plotfun = "stat", stat = "mean")
pp_check(B, plotfun = "dens_overlay")


# Log-linear effect of population density
fC <- n ~ log_pop_dens_scale + IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + offset(log(E))
C <- stan_glm.nb(fC, data = dat_fit, seed = 12345)
# plot(C, plotfun = "trace")
# summary(C)

prior_summary(C)
# pp_check(C, plotfun = "stat", stat = "mean")
# pp_check(C, plotfun = "dens_overlay")


# Saturated effect of population density
## Fit saturation function: Manually optimise saturation function with respect to LOOIC of GLM adjusted for other covariates:

fit_glm <- function(theta){stan_glm.nb(n ~ IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + scale(sat_func(pop_dens, theta)) + offset(log(E)), data = dat_opt, seed = 12345) }

### First set of values
theta <- c(seq(0.001,1,0.001))
fits <- lapply(theta, fit_glm)

loo1 <- lapply(fits, function(fit) loo(fit, save_psis = T)) # retain whole LOO object to check reliability of PSIS
looic1 <- unlist(lapply(loo1, function(loo) loo$estimates[3,1])) # extract metric of interest
png(filename = paste0(figdir, "./with_time_lag/theta_opt_holdout40_totdens.png"), height = 6, width = 8, units = "in", res = 600)
plot(theta, looic1, main = paste0("Optimal theta = ",theta[which.min(looic1)]))
abline(v = theta[which.min(looic1)], col = "red", lty = "dashed")
abline(h = min(looic1), col = "red", lty = "dashed")
dev.off()

# ### Second set
# # theta2 <- c(seq(0.01,0.1,0.01))
# theta2 <- c(seq(0.011,0.059,0.001))
# fits2 <- lapply(theta2, fit_glm)
# looic2 <- unlist(lapply(fits2, function(fit) loo(fit)$estimates[3,1]))
# 
# plot(c(theta, theta2), c(looic1,looic2), main = paste0("Optimal theta = ",theta2[which.min(looic2)]))
# abline(v = theta2[which.min(looic2)], col = "red", lty = "dashed")
# abline(h = min(looic2), col = "red", lty = "dashed")
# 
# # ### Third set
# theta3 <- c(seq(0.001,0.02,0.001))
# # theta3 <- c(seq(0.08,0.099,0.001))
# fits3 <- lapply(theta3, fit_glm)
# looic3 <- sapply(fits3, function(fit) loo(fit)$estimates[3,1])
# # 
# # ### Bind all proposed values and corresponding LOOIC
# theta_all <- c(theta, theta2,theta3)
# fits_all <- rlist::list.append(fits,fits2,fits3)
# looic_all <- unlist(c(looic1,looic2,looic3))
# 
# png(filename = paste0(figdir, "theta_opt_holdout40_totdens_epiwk.png"), height = 6, width = 8, units = "in", res = 600)
# plot(theta_all, looic_all, main = paste0("Optimal theta = ",theta_all[which.min(looic_all)]), xlab = "Theta",ylab = "LOOIC")
# abline(v = theta_all[which.min(looic_all)], col = "red", lty = "dashed")
# abline(h = min(looic_all), col = "red", lty = "dashed")
# dev.off()


## Define theta as that which maximised elpd
# theta_opt <- theta_all[which.min(looic_all)]
theta_opt <- theta[which.min(looic1)]


## Extract optimal fit
# fit_opt <- fits3[[which.in(looic_all)]]
# summary(fit_opt, probs = c(0.01, 0.99), digits = 2)
# 
# prior_summary(fit_opt)
# pp_check(fit_opt, plotfun = "stat", stat = "mean")
# pp_check(fit_opt, plotfun = "dens_overlay")
# posterior_vs_prior(fit_opt)


# Refit final saturating model with main dataset
fD <- n ~ IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + scale(sat_func(pop_dens, theta_opt)) + offset(log(E))
D <- stan_glm.nb(fD, data = dat_fit, seed = 12345)

# Refit without scaling to get coefficient interpretation on original scale
fD2 <- n ~ IMD_scale + prop_minority_scale + prop_kw_scale + wk_since_first_scale + sat_func(pop_dens, theta_opt) + offset(log(E))
D2 <- stan_glm.nb(fD2, data = dat_fit, seed = 12345)
summary(D2)

exp(posterior_interval(D2, pars = c("prop_minority_scale","IMD_scale","prop_kw_scale","wk_since_first_scale","sat_func(pop_dens, theta_opt)"))) #"wk_since_first_scale",
# Without epidemic week:
#                                      5%       95%
# prop_minority_scale           1.0619915 1.2142565
# IMD_scale                     0.9715352 1.0798000
# prop_kw_scale                 0.9272041 1.0137195
# wk_since_first_scale          0.8434189 0.9413276
# sat_func(pop_dens, theta_opt) 1.8755863 3.9506235

# With epidemic week:
#                                      5%       95%
# prop_minority_scale           1.0888155 1.2273346
# IMD_scale                     0.9750827 1.0823019
# prop_kw_scale                 0.9404779 1.0375984
# wk_since_first_scale          0.8575223 0.9577227
# sat_func(pop_dens, theta_opt) 2.1582733 7.0782069

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON
# ---------------------------------------------------------------------------- #

mods <- list(A,B,C,D)
# saveRDS(mods, file = "./mods.rds")
names(mods) <- c("Independent","Linear","Log-linear","Saturating")

## Plot model coefficient posteriors
pdf(file = paste0(figdir,"./with_time_lag/model_coeffs.pdf"), height = 5, width = 6)
lapply(mods, plot_coeffs)
dev.off()

## Calculate approximate LOO CV and compare
loo_mods <- lapply(mods, loo)
comp_mods <- loo_compare(loo_mods) 
print(comp_mods, simplify = F)

looic_mods <- lapply(mods, function(x) loo(x)$estimates[3,])
looictab <- data.frame(Form = names(mods), bind_rows(looic_mods)) %>%
  mutate(Difference = Estimate - min(Estimate)) %>%
  arrange(Difference)

## Plot elpd for all models (maximum is best)
png(filename = paste0(figdir, "./with_time_lag/model_diff_elpd_wwklag.png"), height = 6, width = 8, units = "in", res = 600)
comp_mods %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form") %>%
  filter(Form != "Saturating") %>%
  ggplot(aes(Form, elpd_diff)) + 
  geom_errorbar(aes(ymin = elpd_diff - 2*se_diff, ymax = elpd_diff + 2*se_diff), width = 0.1) +
  geom_point() + 
  theme_classic() + 
  geom_hline(aes(yintercept = 0), lty = "dashed", col = "red") + 
  labs(y = "Difference ELPD", subtitle = "Models compared on expected log predicted density (ELPD), relative to that with the highest value (saturating)\nError bars indicate +/-2SE from the point estimate.")
dev.off()

## Equivalently plot looic for all models (minimum is best)
png(filename = paste0(figdir, "./with_time_lag/model_diff_looic_wwklag.png"), height = 6, width = 8, units = "in", res = 600)
comp_mods %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form") %>%
  filter(Form != "Saturating") %>%
  ggplot(aes(x = Form, y = looic)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = min(looictab$Estimate) - 2*looictab$SE[which.min(looictab$Estimate)], 
                  ymax = min(looictab$Estimate) + 2*looictab$SE[which.min(looictab$Estimate)]), 
              fill= "red",alpha = 0.01) +
  geom_hline(aes(yintercept = min(looictab$Estimate)), lty = "dashed", col = "red") +
  geom_errorbar(aes(ymin = looic - 2*se_looic, ymax = looic + 2*se_looic), width = 0.1) +
  geom_point() + 
  theme_classic() + 
  labs(y = "LOOIC", subtitle = "Red dashed line at LOOIC of preferred (saturating) model. Red rectangle indicates +/- 2SE from this value.\nError bars indicate +/-2SE from the point estimate.")
dev.off()


pred_forms <- dat_fit %>%
  mutate(Saturating = colMeans(posterior_predict(D, offset = log(dat_fit$E))),
         Independent = colMeans(posterior_predict(A, offset = log(dat_fit$E))),
         Linear = colMeans(posterior_predict(B, offset = log(dat_fit$E))),
         Loglinear = colMeans(posterior_predict(C, offset = log(dat_fit$E)))) %>%
  dplyr::select(n, E, pop_dens, Saturating,Independent, Linear, Loglinear) %>%
  pivot_longer(-n:-pop_dens, names_to = "Form") %>%
  arrange(Form)

looictab <- arrange(looictab, Form)
labels <- paste(looictab$Form," (",round(looictab$Estimate,0),")")

png(filename = paste0(figdir, "./with_time_lag/func_forms_predsmr_wwklag.png"), height = 6, width = 8, units = "in", res = 600)
ggplot(pred_forms) +
  geom_smooth(aes(pop_dens,value/E, colour = Form, fill = Form), alpha = 0.1, method = "loess") + 
  geom_point(aes(pop_dens, n/E), alpha = 0.1) +
  geom_smooth(aes(pop_dens, n/E), alpha = 0.1, col = "black", lty = "dashed", method = "loess") +
  scale_x_continuous(trans = "log2") +
  scale_colour_discrete(name = "Form (LOOIC)", labels = labels) +
  # scale_y_continuous(trans = "log") +
  labs(x = "Population density",y = "SMR", 
       title = "Population density versus model-predicted SMR, under four functional forms",
       subtitle = paste0("For saturating function, fitted theta = ",theta_opt,". Black dashed line shows smoother through \nobserved points.")) +
  guides(fill = F) +
  theme_minimal() +
  theme(legend.position = c(0.2,0.8)) 
dev.off()


# ---------------------------------------------------------------------------- #
# IMPACT OF DENSITY REDUCTION
# ---------------------------------------------------------------------------- #

# Refit final saturating model with total dataset
D_tot <- stan_glm.nb(fD, data = dat_tot, seed = 12345)

## Calculate the % difference in predicted deaths according to saturating model, for an 84% reduction in density.
dat_tot_pred <-
  dat_tot %>%
  mutate(red_dens = pop_dens*0.16, # effective density reduction during lockdown
         pred_ini = colMeans(posterior_predict(D_tot, offset = log(dat_tot$E))), 
         pred_red = colMeans(posterior_predict(D_tot, newdata = mutate(dat_tot,pop_dens = red_dens), offset = log(dat_tot$E))),
         perc_red = (pred_ini - pred_red)/pred_ini) 

png(filename = paste0(figdir, "./with_time_lag/saturation_map_wwklag.png"), height = 8, width = 8, units = "in", res = 600)
dat_tot_pred %>%
  full_join(regions, by = "lad19cd") %>%
  basic_map(fill = "perc_red") +
  theme(legend.position = c(0.1,0.8)) +
  labs(title = "Percentage reduction in predicted deaths given \n84% reduction in effective population density ", 
       fill = "")
dev.off()

## Predicted death reduction versus initial population density (by % minority and IMD)
ggplot(dat_fit_pred, aes(pop_dens, perc_red)) + 
  geom_point() + 
  theme_classic() 

mino <- 
  ggplot(dat_fit_pred, aes(pop_dens, perc_red, col = prop_minority)) + 
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c()

imd <- 
  ggplot(dat_fit_pred, aes(pop_dens, perc_red, col = IMD)) + 
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c()

mino + imd


save.image("rerun_totdens_wlag.RData")
