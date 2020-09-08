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


# ---------------------------------------------------------------------------- #
# READ DATA
# ---------------------------------------------------------------------------- #

# Primary dataset:
dat_tot <- readRDS("dat_tot.rds")

## Shapefiles
regions <- readRDS("shp.rds") %>%
  filter(regions, grepl("E", lad19cd)) # filter to England 

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

dat_tot %>%
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = T) +
  labs(title = "COVID-19-associated mortality by LTLA",
       subtitle = paste0("Total deaths per 100,000 between week beginning",min(dat$wod)," and ",max(dat$wod)),
       fill = "") -> map_tot
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

# png(filename = "./figures/map_mort_dens.png", height = 4, width = 8, units = "in", res = 600)
map_SMR + map_dens
# dev.off()

# ---------------------------------------------------------------------------- #
# FITTING
# ---------------------------------------------------------------------------- #

# No effect of population density
fA <- n ~ IMD_scale + prop_minority_scale + prop_kw_scale + offset(log(E))
A <- stan_glm.nb(fA, data = dat_tot, seed = 12345)
# plot(A, plotfun = "trace")
# summary(A)

prior_summary(A)
# pp_check(A, plotfun = "stat", stat = "mean")
# pp_check(A, plotfun = "dens_overlay")


# Linear effect of population density
fB <- n ~ pop_dens_scale + IMD_scale + prop_minority_scale + prop_kw_scale + offset(log(E))
B <- stan_glm.nb(fB, data = dat_tot, seed = 12345)
# plot(B, plotfun = "trace")
# summary(B)

prior_summary(B)
pp_check(B, plotfun = "stat", stat = "mean")
pp_check(B, plotfun = "dens_overlay")


# Log-linear effect of population density
fC <- n ~ logpopdens_scale + IMD_scale + prop_minority_scale + prop_kw_scale + offset(log(E))
C <- stan_glm.nb(fC, data = dat_tot, seed = 12345)
# plot(C, plotfun = "trace")
# summary(C)

prior_summary(C)
# pp_check(C, plotfun = "stat", stat = "mean")
# pp_check(C, plotfun = "dens_overlay")


# Saturated effect of population density
## Fit saturation function: Manually optimise saturation function with respect to LOOIC of GLM adjusted for other covariates:

fit_glm <- function(theta){stan_glm.nb(n ~ IMD_scale + prop_minority_scale + prop_kw_scale + sat_func(pop_dens, theta) + offset(log(E)), data = dat_tot, seed = 12345) }

### First set of values
theta <- c(seq(0.1,1,0.1))
fits <- lapply(theta, fit_glm)

loo1 <- lapply(fits, function(fit) loo(fit, save_psis = T)) # retain whole LOO object to check reliability of PSIS
looic1 <- unlist(lapply(loo1, function(loo) loo$estimates[3,1])) # extract metric of interest
plot(theta, looic1, main = paste0("Optimal theta = ",theta[which.min(looic1)]))

### Second set
theta2 <- c(seq(0.01,0.1,0.01))
fits2 <- lapply(theta2, fit_glm)
looic2 <- lapply(fits2, function(fit) loo(fit)$estimates[3,1])
plot(c(theta,theta2), c(looic1,looic2), main = paste0("Optimal theta = ",theta2[which.min(looic2)]))

### Third set
theta3 <- c(seq(0.025,0.035,0.001))
fits3 <- lapply(theta3, fit_glm)
looic3 <- sapply(fits3, function(fit) loo(fit)$estimates[3,1])

### Bind all proposed values and corresponding elpd
theta_all <- c(theta, theta2,theta3)
fits_all <- rlist::list.append(fits,fits2,fits3)
looic_all <- unlist(c(looic1,looic2,looic3))
waic_all <- lapply(fits_all, function(x) loo::waic(x)$estimates[3,])

# png(filename = "./figures/theta_opt_looic.png",height = 400, width = 600)
plot(theta_all, looic_all, main = paste0("Optimal theta = ",theta_all[which.min(looic_all)]), xlab = "Theta",ylab = "LOOIC")
abline(v = theta_all[which.min(looic_all)], col = "red", lty = "dashed")
abline(h = min(looic_all), col = "red", lty = "dashed")
# dev.off()


## Define theta as that which maximised elpd
theta_opt <- theta_all[which.min(looic_all)]


## Extract optimal fit
D <- fits3[[which.in(looic_all)]]
summary(D, probs = c(0.01, 0.99), digits = 2)

prior_summary(D)
pp_check(D, plotfun = "stat", stat = "mean")
pp_check(D, plotfun = "dens_overlay")

posterior_vs_prior(D)

# fD <- n ~ IMD_scale + prop_minority_scale + prop_kw_scale +  sat_func(pop_dens, theta_opt) + offset(log(E))
# D <- stan_glm.nb(fD, data = dat_tot, seed = 12345)
summary(D)

exp(posterior_interval(D, pars = c("prop_minority_scale","IMD_scale","prop_kw_scale","sat_func(pop_dens, theta)")))

plot(D, pars = c("prop_minority_scale","IMD_scale","prop_kw_scale","sat_func(pop_dens, theta)"))

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON
# ---------------------------------------------------------------------------- #

mods <- list(A,B,C,D)
# saveRDS(mods, file = "./mods.rds")
names(mods) <- c("Independent","Linear","Log-linear","Saturating")

## Plot model coefficient posteriors
# pdf(file = "./figures/model_coeffs.pdf", height = 400, width = 500)
lapply(mods, plot_coeffs)
# dev.off()

## Calculate approximate LOO CV and compare
loo_mods <- lapply(mods, loo)
comp_mods <- loo_compare(loo_mods)
print(comp_mods, simplify = F)

looic_stanmods <- lapply(mods, function(x) loo(x)$estimates[3,])
looictab <- data.frame(Form = names(mods), bind_rows(looic_stanmods)) %>%
  mutate(Difference = Estimate - min(Estimate)) %>%
  arrange(-Difference)

## Plot elpd for all models (maximum is best)
# png(file = "./figures/model_diff.png", height = 400, width = 500)
comp_mods %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form") %>%
  ggplot(aes(Form, elpd_diff)) + 
  geom_errorbar(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), width = 0.1) +
  geom_point() + 
  theme_classic() + 
  geom_hline(aes(yintercept = 0), lty = "dashed", col = "red")
# dev.off()

## Equivalently plot looic for all models (minimum is best)
# png(file = "./figures/model_diff.png", height = 400, width = 500)
looictab %>%
  as.data.frame() %>%
  ggplot(aes(Form, Estimate)) + 
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.1) +
  geom_point() + 
  theme_classic() + 
  geom_hline(aes(yintercept = min(looictab$Estimate)), lty = "dashed", col = "red")
# dev.off()


pred_forms <- dat_tot %>%
  mutate(Saturating = colMeans(posterior_predict(D, offset = log(dat_tot$E))),
         Independent = colMeans(posterior_predict(A, offset = log(dat_tot$E))),
         Linear = colMeans(posterior_predict(B, offset = log(dat_tot$E))),
         Loglinear = colMeans(posterior_predict(C, offset = log(dat_tot$E)))) %>%
  dplyr::select(n, E, pop_dens, Saturating,Independent, Linear, Loglinear) %>%
  pivot_longer(-n:-pop_dens, names_to = "Form")

labels <- paste(sort(unique(dens_forms$Form))," (",round(looictab$Estimate,0),")")

# png(filename = "./figures/func_forms_predsmr.png", height = 6, width = 8, units = "in", res = 600)
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
# dev.off()


# ---------------------------------------------------------------------------- #
# IMPACT OF DENSITY REDUCTION
# ---------------------------------------------------------------------------- #

## Calculate the % difference in predicted deaths according to saturating model, for an 84% reduction in density.
dat_tot_pred <-
  dat_tot %>%
  mutate(red_dens = pop_dens*0.16, # effective density reduction during lockdown
         pred_ini = colMeans(posterior_predict(D, offset = log(dat_tot$E))), 
         pred_red = colMeans(posterior_predict(D, newdata = mutate(dat_tot,pop_dens = red_dens), offset = log(dat_tot$E))),
         perc_red = (pred_ini - pred_red)/pred_ini) 

# png(filename = "./figures/saturation_map.png", height = 8, width = 8, units = "in", res = 600)
dat_tot_pred %>%
  full_join(regions) %>%
  basic_map(fill = "perc_red") +
  theme(legend.position = c(0.1,0.8)) +
  labs(title = "Percentage reduction in predicted deaths given \n84% reduction in effective population density ", 
       fill = "")
# dev.off()

## Predicted death reduction versus initial population density (by % minority and IMD)
ggplot(dat_tot_pred, aes(pop_dens, perc_red)) + 
  geom_point() + 
  theme_classic() 

mino <- 
  ggplot(dat_tot_pred, aes(pop_dens, perc_red, col = prop_minority)) + 
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c()

imd <- 
  ggplot(dat_tot_pred, aes(pop_dens, perc_red, col = IMD)) + 
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c()

mino + imd
