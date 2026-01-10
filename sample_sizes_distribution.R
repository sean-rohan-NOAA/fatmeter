library(akfishcondition)

# Data summary tables

sel_species <- "PCOD"
# sel_species <- "WP"

# Get haul data from RACEBASE
channel <- akfishcondition:::get_connected(schema = "AFSC")

haul <- RODBC::sqlQuery(
  channel = channel,
  query = "SELECT VESSEL, CRUISE, HAUL, HAUL_TYPE, PERFORMANCE, BOTTOM_DEPTH, GEAR_TEMPERATURE, SURFACE_TEMPERATURE, START_LONGITUDE AS LONGITUDE, START_LATITUDE AS LATITUDE FROM RACEBASE.HAUL WHERE REGION = 'BS' and CRUISE > 202100"
)

samples <- readxl::read_xlsx(path = "./data/fatmeter_data_dec_2025.xlsx", sheet = "PCOD") %>%
  dplyr::mutate(`SPECIES NAME` = "Pacific cod",
                `Specimen Number` = as.numeric(`Specimen Number`)) |>
  dplyr::bind_rows(
    readxl::read_xlsx(path = "./data/fatmeter_data_dec_2025.xlsx", sheet = "WP") |>
      dplyr::mutate(`SPECIES NAME` = "walleye pollock")
  ) %>%
  rename_with(~ .x %>%
                str_remove_all("[^[:alnum:] ]") %>%
                str_squish() %>%
                str_replace_all(" ", "_") %>%
                str_to_upper()
  ) |>
  dplyr::mutate(
    YEAR = floor(CRUISE/100),
    PARTIAL_FULLNESS = STOMACH_CONTENT_WT_G/TOTAL_WT_G,
    UNIQUE_HAUL = factor(paste0(VESSEL, "_", CRUISE, "_", HAUL)),
    DOY = lubridate::yday(DATE), # Day of year
    P_LIVLIPID = LIVLIPID/100, # Convert to proportions for model-fitting
    P_MUSLIPID = MUSLIPID/100
  ) |>
  dplyr::select(-LONGITUDE, -LATITUDE) |>
  dplyr::left_join(haul)

table_mean_range <- function(x, digits = 1, nsmall = 1) {
  paste0(format(round(mean(x, na.rm = TRUE), digits = digits), nsmall = nsmall), " (", paste(range(x, na.rm = TRUE), collapse = "-"), ")")
}

sample_table <- 
  samples |>
  dplyr::group_by(SPECIES_NAME, YEAR) |>
  dplyr::summarise(
    n = n(),
    `FM liver` = sum(!is.na(DISTELL_LIVER)),
    `FM muscle` = sum(!is.na(DISTELL_TISSUE)),
    `Ages` = sum(!is.na(AGE)),
    `Liver lipid (%)` = sum(!is.na(LIVLIPID)),
    `Muscle lipid (%)` = sum(!is.na(MUSLIPID)),
    `Ages (yr)` = table_mean_range(AGE),
    `Fork length (cm)` = table_mean_range(LENGTH_CM),
    `Body weight (g)` = table_mean_range(TOTAL_WT_G, digits = 0, nsmall = 0)
  )

map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")


# Exploratory maps ----

# Build summary table for samples and convert to spatial object
spatial_summary <- 
  samples |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, YEAR, LATITUDE, LONGITUDE, SPECIES_NAME) |>
  dplyr::summarise(
    n = n(),
    mean_age = mean(AGE, na.rm = TRUE),
    mean_fl = mean(LENGTH_CM, na.rm = TRUE),
    mean_wt = mean(TOTAL_WT_G, na.rm = TRUE),
    mean_mus_pct_lipid = mean(MUSLIPID, na.rm = TRUE),
    mean_liv_pct_lipid = mean(LIVLIPID, na.rm = TRUE),
    mean_partial_fullness = mean(PARTIAL_FULLNESS, na.rm = TRUE)
  ) |>
  sf::st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338")

p_map_length_pollock <- 
  ggplot() +
  geom_sf(data = map_layers$akland) +
  geom_sf(data = map_layers$survey.area,
          fill = NA, color = "black") +
  geom_sf(
    data = dplyr::filter(spatial_summary, SPECIES_NAME == "walleye pollock"),
    mapping = aes(color = mean_fl, size = n),
    alpha = 0.7
  ) +
  scale_size_continuous(guide = FALSE, range = c(1,4)) +
  scale_x_continuous(
    limits = map_layers$plot.boundary$x,
    breaks = map_layers$lon.breaks
  ) +
  scale_y_continuous(
    limits = map_layers$plot.boundary$y,
    breaks = map_layers$lat.breaks
  ) +
  scale_color_viridis_c(name = "Mean FL (cm)") +
  facet_grid(~YEAR) +
  theme_bw()

p_map_length_pcod <- 
  ggplot() +
  geom_sf(data = map_layers$akland) +
  geom_sf(data = map_layers$survey.area,
          fill = NA, color = "black") +
  geom_sf(
    data = dplyr::filter(spatial_summary, SPECIES_NAME == "Pacific cod"),
    mapping = aes(color = mean_fl, size = n),
    alpha = 0.7
  ) +
  scale_size_continuous(guide = FALSE, range = c(1, 4)) +
  scale_x_continuous(
    limits = map_layers$plot.boundary$x,
    breaks = map_layers$lon.breaks
  ) +
  scale_y_continuous(
    limits = map_layers$plot.boundary$y,
    breaks = map_layers$lat.breaks
  ) +
  scale_color_viridis_c(name = "Mean FL (cm)") +
  facet_grid(~YEAR) +
  theme_bw()

sample_map <- 
  cowplot::plot_grid(
    p_map_length_pollock,
    p_map_length_pcod,
    labels = c("A", "B"),
    nrow = 2
  )

ragg::agg_png(filename = here::here("plots", "sample_maps.png"), height = 80, width = 100, units = "mm")


write.csv(sample_table, file = here::here("plots", "sample_size_table.csv"), row.names = FALSE)
