
# Best P_liver lipid models -----
# Outside of R: selected models for each predictor type based on parsimony (using AIC) and diagnostics.
# Added choices to a spreadsheet

# sel_species <- "PCOD"
sel_species <- "WP"


# Make best model subset table
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
    common_name, category, formula, disp, aic, delta_aic, rmse, mre, mae, r2, mean_bias
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
  x = p_livlipid_table,
  file = here::here("plots", "livlipid_betaglm_table.xlsx"),
  row.names = FALSE
)

# Format full model table for supplement


load(here::here("output", paste0("wp_results.rda")))

wp_results$results_liver_lipid$model_list

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