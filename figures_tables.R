library(akgfmaps)
library(ggpp)


# Best P_liver lipid models -----
# Outside of R: selected models for each predictor type based on parsimony (using AIC) and diagnostics.
# Added choices to a spreadsheet

# Load model outputs
load(here::here("output", paste0("wp_results.rda")))
load(here::here("output", paste0("pcod_results.rda")))

# Best models table -----
model_table <- 
  dplyr::bind_rows(
  read.csv(file = here::here("plots", "WP_liver_model_table.csv")) |>
    dplyr::mutate(name_abbv = "WP"),
  read.csv(file = here::here("plots", "PCOD_liver_model_table.csv")) |>
    dplyr::mutate(name_abbv = "PC")
)

best_p_livlipid_models <- readxl::read_xlsx(path = here::here("output", "best_p_lipid_models.xlsx"))


format_value <-
  function(x, digits) {
    format(round(x, digits = digits), nsmall = digits)
  }

p_livlipid_table <- 
  dplyr::inner_join(
  model_table, best_p_livlipid_models
) |>
  dplyr::select(
    common_name, category, model_name, formula, disp, aic, delta_aic, rmse, mre, mae, r2, mean_bias
  ) |>
  dplyr::mutate(
    aic,
    delta_aic,
    rmse = format_value(rmse*100, 1),
    mre = format_value(mre*100, 1),
    mae = format_value(100*mae, 1),
    r2 = format_value(r2, 2),
    mean_bias = format_value(mean_bias, 3)
  )

xlsx::write.xlsx(
  x = dplyr::select(p_livlipid_table, -model_name),
  file = here::here("plots", "livlipid_betaglm_table.xlsx"),
  row.names = FALSE
)


# Format supplement fit tables for each species ----

format_supplement_tables <- function(common_name, fname) {
  
  clean_fixed_effects <- function(formula_str) {
    form <- as.formula(formula_str)
    
    rhs <- formula(form)[[3]]
    
    term_labels <- attr(terms(form), "term.labels")
    
    fixed_only <- term_labels[!grepl("\\|", term_labels)]
    
    if (length(fixed_only) == 0) {
      return("~1")
    } else {
      return(paste0("~", paste(fixed_only, collapse = " + ")))
    }
  }
  
  output <- read.csv(file = fname) |>
    dplyr::mutate(common_name = common_name) |>
    dplyr::select(
      common_name, model_name, fixed = formula, disp, npar = k, aic, delta_aic, rmse, mre, mae, r2, bias = mean_bias
    ) |>
    dplyr::mutate(
      fixed,
      aic,
      delta_aic,
      npar,
      rmse = format_value(rmse*100, 1),
      mre = format_value(mre*100, 1),
      mae = format_value(100*mae, 1),
      r2 = format_value(r2, 2),
      bias = format_value(bias, 3)
    )
  
  for(ii in 1:nrow(output)) {
    output$fixed[ii] <- clean_fixed_effects(output$fixed[ii])
    
  }
  
  write.csv(x = output, file = gsub(x = fname, pattern = ".csv", replacement = "_formatted.csv"), row.names = FALSE)
  
  return(output)
}

format_supplement_tables(here::here("plots", "WP_liver_model_table.csv"), common_name = "walleye pollock")
format_supplement_tables(here::here("plots", "WP_muscle_model_table.csv"), common_name = "walleye pollock")
format_supplement_tables(here::here("plots", "PCOD_liver_model_table.csv"), common_name = "Pacific cod")
format_supplement_tables(here::here("plots", "PCOD_muscle_model_table.csv"), common_name = "Pacific cod")



# Predictons for the best-fit model ----
# Setup newdata for predictions
make_prediction_vars <- function(x) {
  output <- 
    list(
      Null = 
        expand.grid(
          LENGTH_CM = 
            seq(
              min(x$LENGTH_CM),
              max(x$LENGTH_CM),
              by = 1
            ),
          SPECIES = unique(x$SPECIES),
          UNIQUE_HAUL = factor("dummy"),
          YEAR_FAC = factor("dummy")
        ),
      HSI = 
        expand.grid(
          LENGTH_CM = 
            seq(
              min(x$LENGTH_CM),
              max(x$LENGTH_CM),
              by = 1
            ),
          HSI_PCT = 
            seq(
              min(x$HSI_PCT),
              max(x$HSI_PCT),
              length = 200
            ),
          SPECIES = unique(x$SPECIES),
          UNIQUE_HAUL = factor("dummy"),
          YEAR_FAC = factor("dummy")
        ),
      Morphometric = 
        expand.grid(
          LENGTH_CM = 
            seq(
              min(x$LENGTH_CM),
              max(x$LENGTH_CM),
              by = 1
            ),
          LOG_LW_RESID = 
            seq(
              min(x$LOG_LW_RESID),
              max(x$LOG_LW_RESID),
              length = 200
            ),
          SPECIES = unique(x$SPECIES),
          UNIQUE_HAUL = factor("dummy"),
          YEAR_FAC = factor("dummy")
        ),
      Fatmeter = 
        expand.grid(
          LENGTH_CM = 
            seq(
              min(x$LENGTH_CM),
              max(x$LENGTH_CM),
              by = 1
            ),
          DISTELL_LIVER = 
            seq(
              min(x$DISTELL_LIVER),
              max(x$DISTELL_LIVER),
              length = 200
            ),
          SPECIES = unique(x$SPECIES),
          UNIQUE_HAUL = factor("dummy"),
          YEAR_FAC = factor("dummy")
        )
    )
  
  return(output)
}
best_model_fits <- data.frame()


for(ii in 1:nrow(best_p_livlipid_models)) {
  
  sel_results <- NA
  
  if(best_p_livlipid_models$name_abbv[ii] == "PC") {
    sel_results <- pcod_results
  }
  
  if(best_p_livlipid_models$name_abbv[ii] == "WP") {
    sel_results <- wp_results
  }
    
  pred_vars <- make_prediction_vars(sel_results$dat_complete_liver)
  
  best_model_fits <-
    predict_fits(
      model = sel_results$results_liver_lipid$model_list[[best_p_livlipid_models$model_name[ii]]],
      newdata = pred_vars[[best_p_livlipid_models$category[ii]]],
      # re.form = NULL,
      allow.new.levels = TRUE
    ) |>
    dplyr::mutate(
      common_name = best_p_livlipid_models$common_name[ii],
      model_name = best_p_livlipid_models$model_name[ii],
      category = best_p_livlipid_models$category[ii]
    ) |>
    dplyr::bind_rows(best_model_fits)

  
    
}

# best_model_fits |>
#   dplyr::filter(LENGTH_CM %in% c(50, 60))
# 
# 
# ggplot(data = dplyr::filter(best_model_fits, LENGTH_CM %in% seq(30, 100, 10), category == "HSI")) +
#   geom_ribbon(mapping = aes(x = HSI_PCT, ymin = 100*fit_lwr, ymax = 100*fit_upr), alpha = 0.3) +
#   geom_path(mapping = aes(x = HSI_PCT, y = 100*fit, color = factor(LENGTH_CM), group = LENGTH_CM)) +
#   scale_color_viridis_d(name = "FL (cm)") +
#   facet_wrap(~common_name, scales = "free_x")
# 
# ggplot(data = dplyr::filter(best_model_fits, LENGTH_CM %in% seq(30, 100, 10), category == "Morphometric")) +
#   geom_ribbon(mapping = aes(x = LOG_LW_RESID, ymin = 100*fit_lwr, ymax = 100*fit_upr, fill = LENGTH_CM, group = LENGTH_CM), alpha = 0.3) +
#   geom_path(mapping = aes(x = LOG_LW_RESID, y = 100*fit, color = factor(LENGTH_CM), group = LENGTH_CM)) +
#   scale_color_viridis_d(name = "FL (cm)") +
#   facet_wrap(~common_name, scales = "free_x")
# 
# ggplot(data = dplyr::filter(best_model_fits, LENGTH_CM %in% seq(30, 100, 10), category == "Fatmeter")) +
#   geom_ribbon(mapping = aes(x = DISTELL_LIVER, ymin = 100*fit_lwr, ymax = 100*fit_upr, group = LENGTH_CM), alpha = 0.3) +
#   geom_path(mapping = aes(x = DISTELL_LIVER, y = 100*fit, color = factor(LENGTH_CM))) +
#   scale_color_viridis_d(name = "FL (cm)") +
#   facet_wrap(~common_name, scales = "free_x")
# 
# ggplot(data = ) +
#   geom_ribbon(mapping = aes(x = DISTELL_LIVER, ymin = 100*fit_lwr, ymax = 100*fit_upr, group = LENGTH_CM), alpha = 0.3) +
#   geom_path(mapping = aes(x = DISTELL_LIVER, y = 100*fit, color = factor(LENGTH_CM))) +
#   scale_color_viridis_d(name = "FL (cm)") +
#   facet_wrap(~common_name, scales = "free_x")


fit_55 <- 
  dplyr::filter(best_model_fits, LENGTH_CM == 55) |> 
  dplyr::select(common_name, LENGTH_CM, DISTELL_LIVER, HSI_PCT, LOG_LW_RESID, category, fit, fit_lwr, fit_upr) |>
  tidyr::pivot_longer(cols = c("DISTELL_LIVER", "HSI_PCT", "LOG_LW_RESID")) |>
  dplyr::filter(!is.na(value))


fit_rug <- 
  dplyr::bind_rows(
  pcod_results$dat_complete_liver, 
  dplyr::mutate(wp_results$dat_complete_liver, SPECIMEN_NUMBER = as.character(SPECIMEN_NUMBER))
  ) |>
    dplyr::select(common_name, LOG_LW_RESID, LIVLIPID, HSI_PCT, DISTELL_LIVER, LIVLIPID, YEAR) |>
  tidyr::pivot_longer(cols = c("DISTELL_LIVER", "HSI_PCT", "LOG_LW_RESID")) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::mutate(
    category = ifelse(name == "DISTELL_LIVER", "Fatmeter", ifelse(name == "HSI_PCT", "HSI", "Morphometric"))
  )


p_livlipid_fit <- 
  ggplot() +
  geom_ribbon(data = fit_55, mapping = aes(x = value, ymin = 100*fit_lwr, ymax = 100*fit_upr, group = LENGTH_CM), alpha = 0.3) +
  geom_path(data = fit_55, mapping = aes(x = value, y = 100*fit)) +
  geom_point(
    data = fit_rug,
    mapping = aes(x = value, y = LIVLIPID, fill = factor(YEAR), shape = factor(YEAR)),
    size = 1,
    alpha = 0.7
  ) +
  scale_fill_manual(name = "Year", values = c("#40B0A6", "#5D3A9B")) +
  scale_shape_manual(name = "Year", values = c(21,24)) +
  scale_x_continuous(name = "Predictor value") +
  facet_grid(common_name~category, scales = "free") +
  scale_y_continuous(name = "Liver lipid (%)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(3, units = "mm"),
    legend.key.width = unit(3, units = "mm"),
    legend.spacing = unit(2, units = "mm"),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.12)
  )

png(filename = here::here("plots", "p_livlipid_fit_obs.png"), width = 169, height = 90, units = "mm", res = 300)
print(p_livlipid_fit)
dev.off()
  


hist(residuals(pcod_results$results_liver_lipid$model_list$liver_9))
hist(residuals(wp_results$results_liver_lipid$model_list$liver_18))



# Format full model tables for supplement ----

# TO DO

# Predicted versus observed liver lipid percentage -----


obs_pred <- 
  dplyr::bind_rows(
    wp_results$results_liver_lipid$fits |>
      dplyr::mutate(SPECIMEN_NUMBER = as.character(SPECIMEN_NUMBER)),
    pcod_results$results_liver_lipid$fits
    ) |>
  dplyr::inner_join(best_p_livlipid_models)

obs_pred_labels <-
  p_livlipid_table |>
  dplyr::mutate(
    label = paste0("RMSE=", rmse,"%",
                  "\nMRE=", mre,"%",
                  "\nMAE=", mae,"%",
                  "\nr2=",r2)
  )

p_cond_pred_obs <- 
  ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(
    data = dplyr::filter(obs_pred, category != "Null"),
    mapping = aes(x = LIVLIPID, y = fit*100, fill = factor(YEAR), shape = factor(YEAR)),
    alpha = 0.7
  ) +
  geom_rug(
    data = dplyr::filter(obs_pred, category != "Null"),
    mapping = aes(x = LIVLIPID),
    color = "grey70"
  ) +
  geom_text_npc(
    data = dplyr::filter(obs_pred_labels, category != "Null"),
    mapping = 
      aes(
        npcx = 0.02,
        npcy = 0.98,
        label = label),
    size = 2.2
  ) +
  facet_grid(common_name~factor(paste0(category, " GLM"), levels = c("Null GLM", "HSI GLM", "Fatmeter GLM", "Morphometric GLM"))) +
  scale_y_continuous(name = "Predicted liver lipid (%)") +
  scale_x_continuous(name = "Observed liver lipid (%)") +
  scale_fill_manual(name = "Year", values = c("#40B0A6", "#5D3A9B")) +
  scale_shape_manual(name = "Year", values = c(21,24)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(3, units = "mm"),
    legend.key.width = unit(3, units = "mm"),
    legend.spacing = unit(2, units = "mm"),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.12)
  )


png(filename = here::here("plots", "p_livlipid_obs_vs_pred.png"), width = 169, height = 90, units = "mm", res = 300)
print(p_cond_pred_obs)
dev.off()

# Correlations among different condition metrics ----

# Custom function to calculate r^2 and p-value
plot_cor_grid <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  test <- cor.test(x, y)
  r_sq <- test$estimate^2
  p_val <- test$p.value
  
  sig <- symnum(p_val, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", ".", "ns"))
  
  lbl <- paste0("r²= ", round(r_sq, 2), "\n",
                "p", ifelse(p_val < 0.001, "<0.001", paste0("=", format(p_val, digits = 2))), "\n",
                "(", sig, ")")
  
  ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = lbl, size = 4) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = "grey90"))
}

df <- pcod_results$dat_complete_liver |>
  dplyr::select(
    HSI = HSI_PCT, 
    `Liv. lipid` = LIVLIPID, 
    NLE = RELATIVE_LIVER_ENERGY, 
    TLE = TOTAL_LIVER_ENERGY,
    MCI = LOG_LW_RESID, 
    RCI = RELATIVE_CONDITION
    )

# Generate the plot
p_cor_grid <- 
  ggpairs(df,
        upper = list(continuous = plot_cor_mat),
        lower = list(continuous = wrap("points", alpha = 0.6, color = "steelblue")),
        diag = list(continuous = wrap("densityDiag", fill = "steelblue", color = NA))
) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(3, units = "mm"),
    legend.key.width = unit(3, units = "mm"),
    legend.spacing = unit(2, units = "mm"),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.12)
  )


png(filename = here::here("plots", "indicator_cor_plot.png"), width = 169, height = 169, units = "mm", res = 300)
print(p_cor_grid)
dev.off()

  