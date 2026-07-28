library(akgfmaps)
library(readxl)
library(dplyr)
library(stringr)
library(mgcv)
library(DHARMa)
library(plotly)
library(ggpp)
library(glmmTMB)


# Functions ----------------------------------------------------------------------------------------

# Calculate total liver energy function ----

calc_liver_energy <- function(prop_lipid, liver_mass, lean_protein_ed = 23.6, lipid_ed = 39.5, 
                              include_protein = FALSE) {
    
  lipid_mass <- liver_mass * prop_lipid
  other_mass <- liver_mass - lipid_mass
  lipid_energy <- lipid_mass * lipid_ed
  
  # Option to include protein. In gadids, protein is structural and won't be catabolized for energy
  if(include_protein == TRUE) {
    lean_protein_energy <- other_mass * lean_protein_ed
  } else {
    lean_protein_energy <- 0
  }
  
  total_energy <- lipid_energy + lean_protein_energy
  
  energy_density <- total_energy/liver_mass  

  
  return(
    data.frame(
      TOTAL_LIPID_ENERGY = lipid_energy,
      TOTAL_LIVER_ENERGY = total_energy, 
      LIVER_ENERGY_DENSITY = energy_density
    )
  )
  
}

# Cross validation functions -----------------------------------------------------------------------
run_loocv <- function(model_list, dat) {
  
  # Setup dummy variables for prediction
  original_year <- dat$YEAR_FAC
  original_haul <- dat$UNIQUE_HAUL
  
  dat$UNIQUE_HAUL <- factor("dummy")
  dat$YEAR_FAC <- factor("dummy")
  
  results_list <- lapply(names(model_list), function(m_name) {
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
        re.form = NULL,
        allow.new.levels = TRUE # Allow for new levels
      )
    }
    
    obs_var <- obs <- dat$P_LIVLIPID
    
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
    
    # Replace dummy variables with original
    dat$UNIQUE_HAUL <- dat$UNIQUE_HAUL
    dat$YEAR_FAC <- dat$YEAR_FAC
    
    output <- 
      list(
        perf_metrics = 
          data.frame(
            model_name = m_name, 
            rmse = rmse, 
            mre = mre,
            mae = mae,
            r2 = r2,
            mean_bias = mean_bias
          ),
        fit = 
          cbind(
            dat, 
            data.frame(
              model_name = m_name,
              fit = pred
            )
          )
      )

    
  })
  
  loocv_results <- 
    lapply(
      results_list,
      FUN = function(x) { x$perf_metrics}) |>
    do.call(what = rbind)
  
  oos_fit <-
    lapply(
      results_list,
      FUN = function(x) { x$fit}) |>
    do.call(what = rbind)

  results <-
    list(
      loocv = loocv_results,
      fit = oos_fit
    )
  
  return(results)
}

# Model comparison and diagnostics functions -------------------------------------------------------
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

# Setup data.frames for prediction 

# GLM fitting functions ----------------------------------------------------------------------------
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
  
  # Thin models to ones that passed initial checks
  pass <- aic_table$model_name[aic_table$pass_check == TRUE]
  
  model_list <- model_list[pass]
  
  loocv_out <- run_loocv(model_list = model_list, dat = x)
  
  loocv_table <- loocv_out$loocv
  
  fits <- loocv_out$fit
  
  dharma_results <- dharma_plots(
    model_list = model_list, 
    save_dir = here::here("plots", common_name, paste0(model_name_prefix, "_dharma"))
  )
  
  # Estimate fits
  output <-
    list(
      model_list = model_list,
      aic_table = aic_table,
      fits = fits,
      loocv_table = loocv_table,
      dharma = dharma_results
    )
  
  return(output)
  
} 


# Make predictions ----

predict_fits <- 
  function(model, newdata, re.form = NULL, 
           allow.new.levels = TRUE, ...) {
    
    invlink_fn <- model$modelInfo$family$linkinv
    
    pred <- predict(
      obj = model, 
      newdata = newdata, 
      type = "link",
      re.form = re.form,
      allow.new.levels = allow.new.levels,# Allow for new levels
      se.fit = TRUE
    )
    
    pred <- 
      data.frame(
        fit = invlink_fn(pred$fit),
        fit_lwr = invlink_fn(pred$fit - 2 * pred$se.fit),
        fit_upr = invlink_fn(pred$fit + 2 * pred$se.fit),
        fit_link = pred$fit,
        fit_link_se = pred$se.fit
      )
    
    output <- 
      cbind(
        newdata,
        pred
      )
    
    return(output)
    
  }


# Run analyses -------------------------------------------------------------------------------------

spp_abbv <- c("PCOD", "WP")

for(ii in 1:length(spp_abbv)) {
  
  sel_species <- spp_abbv[ii]
  
  common_name <- sel_common_name <- ifelse(sel_species == "WP", "walleye pollock", "Pacific cod")
  
  sel_species_code <- ifelse(sel_species == "WP", 21740, 21720)
  
  # Set minimum length for calculating LW residuals
  esr_min_length_mm <- 0
  
  # Samples for analysis
  dat_complete_liver <- readRDS(here::here("output", paste0("analysis_samples_", sel_species, ".rds")))
  
  # All samples 
  dat <- readRDS(here::here("output", paste0("all_samples_", sel_species, ".rds")))
  
  p_livlipid_formulas <-
    expand.grid(
      formula = c(
        P_LIVLIPID ~ 1 + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(DISTELL_LIVER, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(DISTELL_LIVER, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(DISTELL_LIVER, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(HSI_PCT, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(HSI_PCT, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(HSI_PCT, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(LOG_LW_RESID, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(LOG_LW_RESID, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_LIVLIPID ~ poly(LOG_LW_RESID, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES)
      ),
      disp = c(~1, ~log(LENGTH_CM))
    )
  
  results_liver_lipid <- 
    fit_beta_glm(
      x = dat_complete_liver, 
      formulas_df = p_livlipid_formulas, 
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
      HSI_PCT = seq(min(dat$HSI_PCT, na.rm = TRUE), max(dat$HSI_PCT, na.rm = TRUE), length = 300),
      LENGTH_CM = min(dat$LENGTH_CM):max(dat$LENGTH_CM), 
      UNIQUE_HAUL = "1", 
      YEAR_FAC = "1"
    )
  
  
  # Fit muscle lipid models --------------------------------------------------------------------------
  dat_complete_muscle <- dplyr::filter(dat_complete_liver, !is.na(DISTELL_TISSUE), !is.na(P_MUSLIPID))
  
  p_muslipid_formulas <-
    expand.grid(
      formula = c(
        P_MUSLIPID ~ 1 + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(DISTELL_TISSUE, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(DISTELL_TISSUE, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(DISTELL_TISSUE, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(LOG_LW_RESID, 1) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(LOG_LW_RESID, 2) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES),
        P_MUSLIPID ~ poly(LOG_LW_RESID, 3) + (1|UNIQUE_HAUL) + (1|YEAR_FAC) + (0 + LENGTH_CM | SPECIES)
      ),
      disp = c(~1, ~log(LENGTH_CM))
    )
  
  results_muscle_lipid <- 
    fit_beta_glm(
      x = dat_complete_liver, 
      formulas_df = p_muslipid_formulas, 
      common_name = common_name,
      model_name_prefix = "muscle"
    )
  
  liver_table <- 
    dplyr::inner_join(
      results_liver_lipid$aic_table,
      results_liver_lipid$loocv_table
    ) |>
    dplyr::filter(pass_check == TRUE)
  
  write.csv(liver_table, file = here::here("plots", paste0(sel_species, "_liver_model_table.csv")), row.names = FALSE)
  
  muscle_table <- 
    dplyr::inner_join(
      results_muscle_lipid$aic_table,
      results_muscle_lipid$loocv_table
    ) |>
    dplyr::filter(pass_check == TRUE)
  
  write.csv(muscle_table, file = here::here("plots", paste0(sel_species, "_muscle_model_table.csv")), row.names = FALSE)
  
  # Save output
  
  assign(
    paste0(tolower(sel_species),"_results"),
    value = list(
      common_name = common_name,
      data = dat,
      dat_complete_liver = dat_complete_liver,
      dat_complete_muscle = dat_complete_muscle,
      results_muscle_lipid = results_muscle_lipid,
      results_liver_lipid = results_liver_lipid
    )
  )
  
  save(
    list = paste0(tolower(sel_species),"_results"),
    file = here::here("output", paste0(tolower(sel_species),"_results.rda"))
  )
  
}

