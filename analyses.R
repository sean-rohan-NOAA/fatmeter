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


# Calculate total liver energy ---------------------------------------------------------------------

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

dat <- cbind(
  dat, 
  calc_liver_energy(prop_lipid = dat$P_LIVLIPID, liver_mass = dat$LIVER_WT_G)
)

# Relative liver energy
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

# Lipid models -------------------------------------------------------------------------------------

# Cross validation ---------------------------------------------------------------------------------
run_loocv <- function(model_list, dat) {
  
  results_list <- lapply(seq_along(model_list), function(m_name) {
    mod <- model_list[[m_name]]
    n_obs <- nrow(dat)
    
    pred <- numeric(n_obs)
    
    for(jj in 1:n_obs) {
      # Update model excluding one observation
      fit_loocv <- update(mod, data = dat[-jj, , drop = FALSE])
      
      # Use mapper function to get the back-transformed prediction
      pred[jj] <- predict(
        obj = fit_loocv, 
        newdata = dat[jj, , drop = FALSE], 
        type = "response",
        re.form = NA, # No random effects-- population level
        allow.new.levels = TRUE # Allow for new levels
      )
    }
    
    obs_var <- 
    
    obs <- dat$P_LIVLIPID
    
    # Root mean square error
    rmse <- sqrt(mean((pred - obs)^2))
    
    # Mean relative error
    mre <- mean((abs(pred - obs))/obs)
    
    # Mean absolute error
    mae <- mean(abs(pred - obs))
    
    # R-squared
    r2 <- cor(pred, obs)^2
    
    # Mean bias
    mean_bias <- mean(pred-obs)
    
    data.frame(
      model_name = m_name, 
      rmse = rmse, 
      mre = mre,
      mae = mae,
      r2 = r2,
      mean_bias = mean_bias
      )
    
  })
  
  results <- do.call(rbind, results_list)
  
  return(results)
}
# Diagnostics --------------------------------------------------------------------------------------
make_aic_table <- 
  function(model_list) {
    
    results <- data.frame(
      model_name = names(model_list) %||% seq_along(model_list),
      formula    = sapply(model_list, function(m) paste(format(formula(m)), collapse = "")),
      
      # Dispersion formula from glmmTMB
      disp       = sapply(
        model_list, 
        function(m){
          if(is(m, "lm")) {
            out <- NA} else{
              out <- paste(format(m$modelInfo$allForm$dispformula))
            }
          out
        }),
      aic        = round(sapply(model_list, AIC), 2),
      k          = sapply(model_list, function(m) attr(logLik(m), "df")),
      convergence = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- m$fit$convergence
            }
          out
        }),
      pdhess = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- m$sdr$pdHess
            }
          out
        }),
      max_gradient = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- max(abs(m$sdr$gradient.fixed))
            }
          out
        }),
      stringsAsFactors = FALSE
    )
    
    results$pass_check <- 
      results$convergence == 0 & results$pdhess & abs(results$max_gradient) < 0.001
    
    results$pass_check[grepl(pattern = "ols", x = results$model_name)] <- TRUE
    
    results$delta_aic <- results$aic - min(results$aic, na.rm = TRUE)
    
    candidates <- results[results$delta_aic < 2, ]
    
    return(results[order(results$aic), ])
    
  } 

dharma_plots <- 
  function(model_list, subset = NULL, nsim = 1000, save_dir) {
    
    model_names    <- names(model_list)
    results        <- vector(mode = "list", length(model_list))
    names(results) <- model_names
    
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    
    for(ii in seq_along(model_list)) {
      
      results[[ii]] <- try(DHARMa::simulateResiduals(model_list[[ii]], n = nsim), silent = TRUE)
      
      if(is(results[[ii]], "try-error")) {
        warning(paste0("warning: Unable to produce DHARMa residuals for model ", model_names[ii]))
        next
      }
      
      fname <- paste0(
        "DHARMa_", gsub(x = model_names[ii], pattern = " ", replacement = "_"), 
        ".png")
      
      png(here::here(save_dir, fname), width = 169, height = 120, units = "mm", res = 300)
      print(plot(results[[ii]]))
      dev.off()
      
    }
    
    return(results)
    
  }

# Fit models ---------------------------------------------------------------------------------------
fit_beta_glm <- function(x, formulas_df, common_name, model_name_prefix = "m") {
  
  model_index <- 1:nrow(formulas_df)
  
  model_list <- 
    lapply(
      model_index, 
      function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = formulas_df$formula[[index]], 
            disp = formulas_df$disp[[index]],
            family = beta_family(link = "logit"), 
            data = x
          )
        
        mod
        
      }
    )
  
  names(model_list) <- paste0(model_name_prefix, "_", model_index)
  
  aic_table <- make_aic_table(model_list)
  
  loocv_table <- run_loocv(model_list = model_list, dat = x)
  
  # Thin models to ones that passed initial checks
  pass <- aic_table$model_name[aic_table$pass_check == TRUE]
  
  model_list <- model_list[pass]
  
  dharma_results <- dharma_plots(
    model_list = model_list, 
    save_dir = here::here("plots", common_name, paste0(model_name_prefix, "_dharma"))
  )
  
  # Estimate fits
  output <-
    list(
      model_list = model_list,
      aic_table = aic_table,
      loocv_table = loocv_table,
      dharma = dharma_results
    )
  
  return(output)
  
} 

# Fatmeter models for liver ----


# Fit liver lipid models ---------------------------------------------------------------------------
dat_complete_liver <- dplyr::filter(dat, !is.na(DISTELL_LIVER), !is.na(P_LIVLIPID))

model_formulas <-
  expand.grid(
    formula = c(
      P_LIVLIPID ~ poly(DISTELL_LIVER, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_LIVLIPID ~ poly(DISTELL_LIVER, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_LIVLIPID ~ poly(DISTELL_LIVER, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_LIVLIPID ~ poly(HSI_PCT, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_LIVLIPID ~ poly(HSI_PCT, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_LIVLIPID ~ poly(HSI_PCT, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES)
    ),
    disp = c(~1, ~log(LENGTH_CM))
  )

results_liver_lipid <- 
  fit_beta_glm(
  x = dat_complete_liver, 
  formulas_df = model_formulas, 
  common_name = common_name,
  model_name_prefix = "liver"
)

# Fork length versus liver lipid for varying levels of Distell fat meter observations
newdata_fatmeter <- 
  expand.grid(
  DISTELL_LIVER = seq(min(dat$DISTELL_LIVER, na.rm = TRUE), max(dat$DISTELL_LIVER, na.rm = TRUE), length = 300),
  LENGTH_CM = min(dat$LENGTH_CM):max(dat$LENGTH_CM), 
  UNIQUE_HAUL = "1", 
  YEAR_FAC = "1"
)

newdata_hsi <- 
  expand.grid(
    DISTELL_LIVER = seq(min(HSI_PCT, na.rm = TRUE), max(HSI_PCT, na.rm = TRUE), length = 300),
    LENGTH_CM = min(dat$LENGTH_CM):max(dat$LENGTH_CM), 
    UNIQUE_HAUL = "1", 
    YEAR_FAC = "1"
  )


# Fit muscle lipid models --------------------------------------------------------------------------
dat_complete_muscle <- dplyr::filter(dat, !is.na(DISTELL_TISSUE), !is.na(P_MUSLIPID))

model_formulas <-
  expand.grid(
    formula = c(
      P_MUSLIPID ~ 1 + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_MUSLIPID ~ poly(DISTELL_TISSUE, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_MUSLIPID ~ poly(DISTELL_TISSUE, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
      P_MUSLIPID ~ poly(DISTELL_TISSUE, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES)
    ),
    disp = c(~1, ~log(LENGTH_CM))
  )

results_muscle_lipid <- 
  fit_beta_glm(
    x = dat_complete_liver, 
    formulas_df = model_formulas, 
    common_name = common_name,
    model_name_prefix = "muscle"
  )

# Save output

assign(
  paste0(tolower(sel_species),"_results"),
  value = list(
    common_name = common_name,
    data = dat,
    results_muscle_lipid = results_muscle_lipid,
    results_liver_lipid = results_liver_lipid,
    lw_model = lw_mod
  )
)

save(
  list = paste0(tolower(sel_species),"_results"),
  file = here::here("output", paste0(tolower(sel_species),"_results.rda"))
  )





p_fl_vs_liver_lipid <- 
  ggplot() +
  geom_path(
    data = dplyr::filter(best_fit_liver, DISTELL_LIVER %%10 == 0),
    mapping = aes(x = LENGTH_CM, y = fit*100, group = DISTELL_LIVER, color = DISTELL_LIVER)
  ) +
  geom_point(
    data = dat,
    mapping = aes(
      x = LENGTH_CM, y = P_LIVLIPID*100, color = DISTELL_LIVER)
  ) +
  scale_fill_viridis_c(name = "Fatmeter", na.value = NA, limits = range(dat$DISTELL_LIVER, na.rm = TRUE)) +
  scale_color_viridis_c(name = "Fatmeter", na.value = NA, limits = range(dat$DISTELL_LIVER, na.rm = TRUE)) +
  scale_x_continuous(name = "Fork length (cm)") +
  scale_y_continuous(name = "Liver lipid (%)") +
  theme_bw()

p_fatmeter_vs_lipid_liver <- 
  ggplot() +
  geom_point(
    data = dplyr::filter(dat, !is.na(LENGTH_CM)),
    mapping = aes(
      x = DISTELL_LIVER, y = P_LIVLIPID*100, color = cut(LENGTH_CM, breaks = seq(0, 100, 10))
    )
  ) +
  geom_path(
    data = dplyr::filter(best_fit_liver, LENGTH_CM %%10 == 5),
    mapping = aes(x = DISTELL_LIVER, y = fit*100, color = cut(LENGTH_CM, breaks = seq(0, 100, 10)))
  ) +
  geom_point(
    data = dplyr::filter(dat, !is.na(LENGTH_CM)),
    mapping = aes(
      x = DISTELL_LIVER, y = P_LIVLIPID*100, color = cut(LENGTH_CM, breaks = seq(0, 100, 10))
    )
  ) +
  scale_color_viridis_d(name = "FL (cm)", na.value = NA, na.translate = FALSE, direction = -1, option = "A") +
  scale_fill_viridis_d(name = "FL (cm)", na.value = NA, na.translate = FALSE, direction = -1, option = "A") +
  scale_x_continuous(name = "Fatmeter") +
  scale_y_continuous(name = "Liver lipid (%)") +
  theme_bw()

p_fl_vs_muscle_lipid <- 
  ggplot() +
  geom_path(
    data = dplyr::filter(best_fit_muscle, DISTELL_TISSUE %%10 == 0),
    mapping = aes(x = LENGTH_CM, y = fit*100, group = DISTELL_TISSUE, color = DISTELL_TISSUE)
  ) +
  geom_point(
    data = dat,
    mapping = aes(
      x = LENGTH_CM, y = P_MUSLIPID*100, color = DISTELL_TISSUE)
  ) +
  scale_x_continuous(name = "Fork length (cm)") +
  scale_y_continuous(name = "Liver lipid (%)") +
  scale_color_viridis_c(name = "Fatmeter", na.value = NA) +
  theme_bw()






