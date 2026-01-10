# Calculate morphometric condition ----

# Retrieve haul data - only need to do this once

# channel <- akfishcondition:::get_connected(schema = "AFSC")
# 
# haul_dat <- 
#   RODBC::sqlQuery(
#   channel = channel,
#   query = "select * from racebase.haul where region = 'BS' and cruise >= 202101"
# ) |>
#   dplyr::mutate(area_id = ifelse(STRATUM < 70, floor(STRATUM / 10), STRATUM))
# 
# saveRDS(object = haul_dat, file = here::here("data", "haul_data.rds"))


# Setup analysis parameters ----




