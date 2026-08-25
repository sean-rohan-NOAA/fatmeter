# Make sample map

library(akgfmaps) # GitHub: afsc-gap-products/akgfmaps
library(ggthemes)
library(shadowtext)

format_value <-
  function(x, digits) {
    format(round(x, digits = digits), nsmall = digits, big.mark = ",", trim = TRUE)
  }

set_table_val <- function(x_mean, x_min, x_max, digits) {
  
  lab <- paste0(
    format_value(x_mean, digits), " (",
    format_value(x_min, digits), "-",
    format_value(x_max, digits), ")"
  )
  
}

dat <- 
  dplyr::bind_rows(
    readRDS(here::here("output", "analysis_samples_WP.rds")) |>
      dplyr::mutate(SPECIMEN_NUMBER = as.character(SPECIMEN_NUMBER)),
    readRDS(here::here("output", "analysis_samples_PCOD.rds"))
  )

# Load haul data and join with fatmeter data

haul_dat <- readRDS(here::here("data", "haul_data.rds")) |>
  dplyr::select(VESSEL, CRUISE, HAUL, area_id, X = START_LONGITUDE, Y = START_LATITUDE)

map_dat <- 
  dplyr::left_join(
    dat,
    haul_dat
  ) |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, YEAR, X, Y, common_name) |>
  dplyr::summarise(n = n()) |> 
  sf::st_as_sf(coords = c("X", "Y"), crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338")

map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

land_label <-
  data.frame(
    label = "Alaska",
    x = -160,
    y = 62
  ) |>
  sf::st_as_sf(
    coords = c("x", "y"),
    crs = "WGS84"
  ) |>
  sf::st_transform(crs = "EPSG:3338")

land_label[c("x", "y")] <- sf::st_coordinates(land_label)

p_sample_map <- 
  ggplot() +
  geom_sf(
    data = map_layers$akland,
    color = NA,
    fill = "grey80"
    ) +
  geom_sf(
    data = map_layers$survey.area,
    fill = NA
  ) +
  geom_sf(
    data = map_layers$bathymetry,
    mapping = aes(color = factor(DEPTH_M)),
    linewidth = 0.4
  ) +
  geom_sf(
    data = map_dat,
    mapping = aes(fill = factor(YEAR),
                  shape = factor(YEAR),
                  size = n),
    alpha = 0.5
  ) +
  geom_sf_text(
    data = land_label,
    mapping = aes(x = x, y = y, label = label),
    color = "black",
    size = 3
  ) +
  scale_shape_manual(name = "Year", values = c(21, 24)) +
  scale_size(name = "Fish (#)", range = c(1,3)) +
  scale_x_continuous(
    limits = map_layers$plot.boundary$x,
    breaks = map_layers$lon.breaks
    ) +
  scale_y_continuous(
    limits = map_layers$plot.boundary$y,
    breaks = map_layers$lat.breaks
    ) +
  scale_fill_manual(name = "Year", values = c("#40B0A6", "#5D3A9B")) +
  scale_color_brewer(name = "Depth (m)") +
  facet_wrap(~common_name) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(3, units = "mm"),
    legend.key.width = unit(3, units = "mm"),
    legend.spacing = unit(2, units = "mm")
        )

png(filename = here::here("plots", "sample_map.png"), width = 140, height = 70, units = "mm",
    res = 300)
print(p_sample_map)
dev.off()


# Summary table -----
summary_dat <- dat |>
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
  dplyr::left_join(
    haul_dat
  )


summary_table <- 
  summary_dat |>
  dplyr::group_by(
    YEAR, common_name
  ) |>
  dplyr::summarise(
    n = n(),
    FL_mean = mean(LENGTH_CM, na.rm = TRUE),
    FL_min = min(LENGTH_CM, na.rm = TRUE),
    FL_max = max(LENGTH_CM, na.rm = TRUE),
    BW_mean = mean(TOTAL_WT_G/1000, na.rm = TRUE),
    BW_min = min(TOTAL_WT_G/1000, na.rm = TRUE),
    BW_max = max(TOTAL_WT_G/1000, na.rm = TRUE),
    AGE_mean = mean(AGE, na.rm = TRUE),
    AGE_min = min(AGE, na.rm = TRUE),
    AGE_max = max(AGE, na.rm = TRUE),
    LIV_mean = mean(LIVER_WT_G),
    LIV_min = min(LIVER_WT_G),
    LIV_max = max(LIVER_WT_G)
    ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FL =  set_table_val(FL_mean, FL_min, FL_max, 1),
    BWT_KG =  set_table_val(BW_mean, BW_min, BW_max, 2),
    AGE =  set_table_val(AGE_mean, AGE_min, AGE_max, 1),
    LIV_G = set_table_val(LIV_mean, LIV_min, LIV_max, 1)
  ) |>
  dplyr::select(common_name, YEAR, n, FL, BWT_KG, AGE, LIV_G) |>
  dplyr::arrange(common_name, YEAR)

xlsx::write.xlsx(
  as.data.frame(summary_table),
  file = here::here("plots", "sample_table.xlsx"),
  row.names = FALSE
)
