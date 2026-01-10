library(akgfmaps)
library(readxl)
library(dplyr)
library(stringr)
library(mgcv)
library(DHARMa)
library(plotly)
library(ggpp)
library(glmmTMB)

sel_species <- "PCOD"
# sel_species <- "WP"

common_name <- ifelse(sel_species == "WP", "walleye pollock", "Pacific cod")

esr_min_length_mm <- ifelse(sel_species == "WP", 250, 0)


# Functions to calculate liver energy density and total energy ----

calc_liver_energy <- function(prop_lipid, liver_mass, lean_protein_ed = 23.6, lipid_ed = 39.5) {
  
  lipid_mass <- liver_mass * prop_lipid
  other_mass <- liver_mass - lipid_mass
  lipid_energy <- lipid_mass * lipid_ed
  lean_protein_energy <- other_mass * lean_protein_ed
  
  total_energy <- lipid_energy + lean_protein_energy
  
  energy_density <- total_energy/liver_mass
  
  return(data.frame(
    TOTAL_LIPID_ENERGY = lipid_energy,
    TOTAL_LIVER_ENERGY = total_energy, 
    LIVER_ENERGY_DENSITY = energy_density
  ))
}

dat <- readxl::read_xlsx(path = "./data/fatmeter_data_dec_2025.xlsx", sheet = sel_species) %>%
  rename_with(~ .x %>%
                str_remove_all("[^[:alnum:] ]") %>%
                str_squish() %>%
                str_replace_all(" ", "_") %>%
                str_to_upper()
  ) |>
  dplyr::mutate(
    LENGTH_MM = LENGTH_CM*10,
    common_name = common_name,
    YEAR = floor(CRUISE/100),
    YEAR_FAC = factor(YEAR),
    PARTIAL_FULLNESS = STOMACH_CONTENT_WT_G/TOTAL_WT_G,
    UNIQUE_HAUL = factor(paste0(VESSEL, "_", CRUISE, "_", HAUL)),
    DOY = lubridate::yday(DATE), # Day of year
    P_LIVLIPID = LIVLIPID/100, # Convert to proportions for model-fitting
    P_MUSLIPID = MUSLIPID/100,
    HSI_PCT = LIVER_WT_G/(TOTAL_WT_G-STOMACH_CONTENT_WT_G-OVARY_WT_G)*100,
    GSI_PCT = OVARY_WT_G/(TOTAL_WT_G-STOMACH_CONTENT_WT_G-LIVER_WT_G)*100
  )


# Calculate total liver energy, liver energy density, total lipid energy, and relative liver energy
dat <- cbind(
  dat, 
  calc_liver_energy(prop_lipid = dat$P_LIVLIPID, liver_mass = dat$LIVER_WT_G)
)

# Relative liver energy
dat$RELATIVE_LIVER_ENERGY <- dat$TOTAL_LIVER_ENERGY/dat$TOTAL_WT_G
dat$RELATIVE_LIPID_ENERGY <- dat$TOTAL_LIPID_ENERGY/dat$TOTAL_WT_G

ggplot() +
  geom_point(
    data = dat,
    mapping = aes(x = HSI_PCT, y = RELATIVE_LIVER_ENERGY)
  )

ggplot() +
  geom_point(
    data = dat,
    mapping = aes(x = HSI_PCT, y = RELATIVE_LIPID_ENERGY)
  )



# Calculate morphometric condition -----------------------------------------------------------------

# Load haul data and join with fatmeter dats

haul_dat <- readRDS(here::here("data", "haul_data.rds")) |>
  dplyr::select(vessel = VESSEL, cruise = CRUISE, haul = HAUL, area_id)

morph_dat <- 
  dat |>
  dplyr::select(
    vessel = VESSEL, 
    cruise = CRUISE, 
    haul = HAUL, 
    common_name,
    year = YEAR, 
    specimenid = SPECIMEN_NUMBER, 
    length_mm = LENGTH_MM, 
    weight_g = TOTAL_WT_G
  ) |>
  dplyr::left_join(haul_dat) |>
  dplyr::mutate(area_id = factor(area_id))

# Get all of the akfishcondition data from sampling years that aren't in fatmeter samples
akfishcondition_dat <- 
  dplyr::bind_rows(
    read.csv(file = here::here("data", "ebs_all_species.csv")),
    read.csv(file = here::here("data", "nbs_all_species.csv"))
  ) |>
  dplyr::filter(sex == 2, year %in% unique(morph_dat$year), common_name == common_name, length_mm >= esr_min_length_mm) |>
  dplyr::mutate(
    length_cm = length_mm/10,
    area_id = factor(area_id),
    specimenid = as.character(specimenid)
  ) |>
  dplyr::select(vessel, cruise, haul, area_id, stratum, year, common_name, specimenid, length_mm, weight_g) |>
  dplyr::anti_join(
    morph_dat, by = c("vessel", "cruise", "haul", "specimenid")
  )

combined_dat <- dplyr::bind_rows(morph_dat, akfishcondition_dat)

# Linear regression between log-length and log-weight


lw_mod <- lm(formula = log(weight_g)~log(length_mm), data = combined_dat)
lw_fit <- predict(lw_mod, newdata = morph_dat, se.fit = TRUE)
sigma2 <- summary(lw_mod)$sigma^2
morph_dat$fit <- exp(lw_fit$fit + 0.5*sigma2)
morph_dat$fit_lwr <- exp(lw_fit$fit - 2*lw_fit$se.fit + 0.5*sigma2)
morph_dat$fit_upr <- exp(lw_fit$fit + 2*lw_fit$se.fit + 0.5*sigma2)

dat$LOG_LW_RESID <- log(morph_dat$weight_g)-log(morph_dat$fit)
dat$RELATIVE_CONDITION <- morph_dat$weight_g/morph_dat$fit

# Exploratory plots
# ggplot(data = dplyr::filter(dat, abs(LOG_LW_RESID) < 0.35),
#        mapping = aes(x = LOG_LW_RESID, y = HSI_PCT)) +
#   geom_point() +
#   geom_smooth(method = 'lm')
# 
# 
# cowplot::plot_grid(
#   ggplot(data = dat,
#          mapping = aes(x = AGE, LIVLIPID)) +
#     geom_point(
#     ) +
#     geom_smooth(method = 'lm'),
#   ggplot(data = dat,
#          mapping = aes(x = TOTAL_WT_G, LIVLIPID)) +
#     geom_point(
#     ) +
#     geom_smooth(method = 'lm')
# )
# 
# 
# ggplot(data = dat,
#        mapping = aes(x = LIVER_WT_G, y = LIVLIPID, color = factor(AGE))) +
#   geom_point() +
#   # geom_smooth() +
#   scale_color_viridis_d(name = "Age") +
#   scale_x_continuous(name = "Liver weight (g)") +
#   scale_y_continuous(name = "Liver lipid (%)")
# 
# ggplot(data = dat,
#        mapping = aes(x = HSI_PCT, y = LIVLIPID, color = factor(AGE))) +
#   geom_point() +
#   # geom_smooth() +
#   scale_color_viridis_d(name = "Age") +
#   scale_x_continuous(name = "HSI (%)") +
#   scale_y_continuous(name = "Liver lipid (%)") +
#   theme_bw()
# 
# ggplot(data = dat,
#        mapping = aes(x = AGE, y = LIVLIPID)) +
#   geom_point() +
#   geom_smooth() +
#   ggtitle(label = sel_species) +
#   theme_bw()
# 
# ggplot(data = dat,
#        mapping = aes(x = LENGTH_CM, y = LIVLIPID)) +
#   geom_smooth() +
#   geom_point(mapping = aes(color = factor(AGE))) +
#   scale_color_viridis_d(name = "Age", na.value = "grey") +
#   scale_x_continuous(name = "Fork length (cm)") +
#   scale_y_continuous(name = "Liver lipid (%)") +
#   ggtitle(label = sel_species) +
#   theme_bw()
# 
# ggplot(data = as.data.frame(dat), mapping = aes(x = LENGTH_CM, y = TOTAL_WT_G, color = factor(AGE))) +
#   geom_point(
#   ) +
#   scale_x_continuous(name = "Fork length (cm)") +
#   scale_y_continuous(name = "Total weight (g)") +
#   scale_color_viridis_d(name = "Age", na.value = "grey") +
#   theme_bw()


# Checking for outliers
plot_ly(
  x = dat$LENGTH_CM, y = dat$TOTAL_WT_G, text =  paste0(dat$CRUISE, " ", dat$VESSEL, " ", dat$HAUL, " ", dat$SPECIMEN_NUMBER)
)


# GAM fit to total weight ----

m0_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM) + s(UNIQUE_HAUL, bs = "re"),
  data = dat,
  family = tw()
)

m1_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM) + s(LIVER_WT_G) + s(UNIQUE_HAUL, bs = "re"),
  data = dat,
  family = tw()
)

m2_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM) + s(STOMACH_CONTENT_WT_G) + s(UNIQUE_HAUL, bs = "re"),
  data = dat,
  family = tw()
)

m3_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM) + s(OVARY_WT_G) + s(UNIQUE_HAUL, bs = "re"),
  data = dat,
  family = tw()
)

m4_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM) + s(LIVER_WT_G) + s(OVARY_WT_G) + s(UNIQUE_HAUL, bs = "re"),
  data = dat,
  family = tw()
)

m5_total_weight <- mgcv::gam(
  formula = TOTAL_WT_G ~ s(LENGTH_CM, bs = "tp") + s(LIVER_WT_G, bs = "tp") + s(OVARY_WT_G, bs = "re") + s(STOMACH_CONTENT_WT_G, bs = "re") + s(UNIQUE_HAUL, bs = "re"),
  family = tw(link = "log"),
  data = dat)

summary(m0_total_weight)
summary(m1_total_weight)
summary(m2_total_weight)
summary(m3_total_weight)
summary(m4_total_weight)
summary(m5_total_weight)

gam.check(m5_total_weight)

aic_total_weight <-
  AIC(
    m0_total_weight,
    m1_total_weight,
    m2_total_weight,
    m3_total_weight,
    m4_total_weight,
    m5_total_weight
  )

# Residuals for the best model
dharma_resids <- DHARMa::simulateResiduals(m5_total_weight, n = 250)