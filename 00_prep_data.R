library(akgfmaps) # GitHub: afsc-gap-products/akgfmaps
library(readxl)
library(dplyr)
library(stringr)
library(mgcv)
library(DHARMa)
library(plotly)
library(ggpp)
library(glmmTMB)

spp_abbv <- c("PCOD", "WP")

for(ii in 1:length(spp_abbv)) {
  
  sel_species <- spp_abbv[ii]
  
  common_name <- sel_common_name <- ifelse(sel_species == "WP", "walleye pollock", "Pacific cod")
  
  sel_species_code <- ifelse(sel_species == "WP", 21740, 21720)
  
  # Set minimum length for calculating LW residuals
  esr_min_length_mm <- 0
  
  ftnir_dat <- readxl::read_xlsx("./data/Prohaska_EFH_liver_NIRpreds.xlsx", sheet = 1) |>
    dplyr::filter(species_code == sel_species_code, region == 'BS') |>
    dplyr::mutate(common_name = sel_common_name)
  
  ftnir_samples <- 
    ftnir_dat |>
    dplyr::select(VESSEL = vessel_code, CRUISE = cruise_number, HAUL = haul, LENGTH_MM = length, TOTAL_WT_G = weight, Lipid_pred_LOOCV, Lipid_pred_cal)
  
  dat <- readxl::read_xlsx(path = "./data/fatmeter_data_feb_2026.xlsx", sheet = sel_species) %>%
    dplyr::filter(!(comments == "exclude") | is.na(comments)) %>%
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
    ) |>
    dplyr::filter(!is.na(LIVLIPID))
  
  # Only use the subset of the data that have corresponding FT-NIR samples
  dat <- dplyr::inner_join(dat, ftnir_samples)
  
  # Check for completeness
  dplyr::anti_join(dat, ftnir_samples)
  dplyr::anti_join(ftnir_samples, dat)
  
  p_lw <- ggplot() +
    geom_point(
      data = dat, 
      mapping = aes(x = LENGTH_CM, y = TOTAL_WT_G), size = 0.2) +
    scale_x_continuous(name = "Fork length (cm)") +
    scale_y_continuous(name = "Total weight (g)") +
    theme_bw()
  
  png(filename = here::here("plots", paste0(sel_species, "_length_weight.png")), width = 80, height = 80, units = "mm", res = 300)
  print(p_lw)
  dev.off()
  
  # Relative liver energy
  dat <- cbind(
    dat, 
    calc_liver_energy(prop_lipid = dat$P_LIVLIPID, liver_mass = dat$LIVER_WT_G)
  )
  
  dat$RELATIVE_LIVER_ENERGY <- dat$TOTAL_LIVER_ENERGY/dat$TOTAL_WT_G
  dat$RELATIVE_LIPID_ENERGY <- dat$TOTAL_LIPID_ENERGY/dat$TOTAL_WT_G
  
  
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
    dplyr::mutate(area_id = factor(area_id),
                  specimenid = as.character(specimenid))
  
  # Get all of the akfishcondition data from sampling years that aren't in fatmeter samples
  akfishcondition_dat <- 
    dplyr::bind_rows(
      read.csv(file = here::here("data", "ebs_all_species.csv")),
      read.csv(file = here::here("data", "nbs_all_species.csv"))
    ) |>
    dplyr::filter(
      sex == 2, 
      year %in% unique(morph_dat$year), 
      common_name == sel_common_name, 
      length_mm >= esr_min_length_mm
    ) |>
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
  
  
  # Fit liver lipid models ---------------------------------------------------------------------------
  dat_complete_liver <- dplyr::filter(dat, !is.na(DISTELL_LIVER), !is.na(P_LIVLIPID))
  
  saveRDS(dat_complete_liver, here::here("output", paste0("analysis_samples_", spp_abbv[ii], ".rds")))
  
  saveRDS(dat, here::here("output", paste0("all_samples_", spp_abbv[ii], ".rds")))
  
}
