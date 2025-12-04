# =============================================================================
# Script name: run_transfer_simulations.R
# Author: Jack Litle
# Date: 2025-12-02
# Purpose: Run biophysical simulations for PNWHD mussel heatwave study
# Notes: Only core simulation pipeline included; figure/statistical scripts not included
# =============================================================================

#Load required packages
library(tidyverse)
library(TrenchR)


# Source function definitions
source("../R/heatbudget_functions.R")

# =========================================================================
# SECTION 1: Load and/or create template data with weather conditions 
# and mussel sizes
# =========================================================================

#load example input data
exposure_df_hourly <- read.csv("../data/example_inputs.csv", header = TRUE) 
head(exposure_df_hourly)

#Example file only includes only one length, expand to many lengths if need be:
#Step 1: Your sequence of lengths
new_lengths <- seq(0.02, 0.200, by = 0.01)

# Step 2: Filter to a single template size (just for environmental conditions)
base_df_single_size <- exposure_df_hourly %>% 
  filter(length == 0.08) %>%        # pick any one size
  select(-length)                   # remove the 'length' column entirely

# Step 3: Repeat rows for each new length
exposure_df_hourly <- base_df_single_size %>% 
  crossing(length = new_lengths)

exposure_df_hourly <- as.data.frame(exposure_df_hourly) #make it into a df 


# ===========================
# SECTION 2: Run Simulations
# ===========================

# Group by length
group_var <- "length"

# Sort by group and aerial exposure time
exposure_df_hourly <- exposure_df_hourly[order(exposure_df_hourly[[group_var]], exposure_df_hourly$exposure_time), ]

# Pre-create necessary columns if they don't exist
temp_cols <- c("downscaled_wind", "solar_zenith_angle", "prop_diffuse_solar",
               "remaining_water", "percent_water_lost", "evap",
               "tissue_mass", "Q_stored", "body_temp", "delta_Tb")

for (col in temp_cols) {
  if (!(col %in% names(exposure_df_hourly))) {
    exposure_df_hourly[[col]] <- NA_real_
  }
}

# Define constants
c_shell <- 815     # J/kg·K
c_body <- 4180     # J/kg·K
T_init <- 30 + 273.15  # K (initial body temp - mean of all the robomussel temps after a time step of aerial exposure)
fade_start <- 0.245 #simulates mussels closing when they are close to water threshold
fade_end <- 0.255 # this is about a 10 percent water loss limit
#fade_start <- 0.535 #20% water loss limit
#fade_end <- 0.54
#fade_start <- 0.805  #30% water loss limit
#fade_end <- 0.810


for (m in unique(exposure_df_hourly[[group_var]])) { #define groups by mussel length
  group_rows <- which(exposure_df_hourly[[group_var]] == m)
  
  #define length, initial water, tissue and shell mass, and initial body temp
  L <- m 
  water_initial <- initial_water(L)
  tissue_fixed <- tissue_mass(L)
  m_shell <- shell_mass(L)
  
  water_remaining <- water_initial
  m_body <- tissue_fixed
  
  T_b <- T_init
  
  cat(sprintf("\n=== Starting group for L = %.5f ===\n", L))
  cat(sprintf("Initial water mass = %.6f kg, Tissue mass = %.6f kg, Shell mass = %.6f kg\n", 
              water_initial, tissue_fixed, m_shell))
  
  for (j in seq_along(group_rows)) { #within each group
    i <- group_rows[j]
    
    #downscale wind, set solar zenith angle, and determine diffuse solar radiation
    exposure_df_hourly[i, "downscaled_wind"] <- wind_speed_profile_neutral(
      u_r = exposure_df_hourly[i, "wind_speed"],
      zr = 2.5, z0 = 0.3, z = 0.5
    )
    
    exposure_df_hourly[i, "solar_zenith_angle"] <- zenith_angle(
      hour = exposure_df_hourly[i, "hour"],
      doy = exposure_df_hourly[i, "doy"],
      lat = 48.4504,
      lon = -122.9635
    )
    
    exposure_df_hourly[i, "prop_diffuse_solar"] <- proportion_diffuse_solar_radiation(
      psi = exposure_df_hourly[i, "solar_zenith_angle"],
      p_a = 101.3, rho = 0.15
    )
    
    #use downscaled wind
    wind <- exposure_df_hourly[i, "downscaled_wind"]
    
    #convert temps to Kelvin (input in Celsius) annd define relative humidity
    T_air <- exposure_df_hourly[i, "air_temp"] + 273.15
    T_g <- exposure_df_hourly[i, "ground_temp"] + 273.15
    RH <- exposure_df_hourly[i, "RH"]
    
    #define vapor densities of body and air
    vd_body <- SatVapDense(T_b)
    vd_air <- SatVapDense(T_air) * RH
    vd_diff <- vd_body - vd_air
    
    #calculate evaporation rate
    evap_raw <- evap_rate(wind = wind, L = L, vd_body = vd_body, vd_air = vd_air)
    if (evap_raw < 0) evap_raw <- 0
    
    # Time step calculation
    if (j == 1) {
      time_step <- exposure_df_hourly[i, "exposure_time"]
    } else {
      t_now <- exposure_df_hourly[i, "exposure_time"]
      t_prev <- exposure_df_hourly[group_rows[j - 1], "exposure_time"]
      time_step <- t_now - t_prev
    }
    
    #calculate cumulative waer loss
    cumulative_water_loss <- water_initial - water_remaining
    p_lost <- cumulative_water_loss / water_initial
    
    # Steep linear fade simulating mussels closing valves when they get to threshold water
    fade_factor <- case_when(
      p_lost <= fade_start ~ 1,
      p_lost >= fade_end ~ 0,
      TRUE ~ (fade_end - p_lost) / (fade_end - fade_start)
    )
    
    evap <- evap_raw * fade_factor
    water_loss <- evap * time_step
    water_remaining <- max(0, water_remaining - water_loss)
    
    m_body <- tissue_fixed - (water_initial - water_remaining)
    percent_lost <- 100 * (water_initial - water_remaining) / water_initial
    
    exposure_df_hourly[i, "remaining_water"] <- water_remaining
    exposure_df_hourly[i, "percent_water_lost"] <- percent_lost
    exposure_df_hourly[i, "evap"] <- fade_factor > 0.001  # logical TRUE/FALSE
    exposure_df_hourly[i, "tissue_mass"] <- m_body
    
    # Define full heat capacity of the system
    C_total <- m_shell * c_shell + c_body * m_body
    
    # Run full heat transfer simulation for timestep
    if (fade_factor > 0.001) { # if the mussel has lost less than 25% of available water
      Q_rate <- Q_unsteady_evap( # use the evap function
        T_a = T_air, T_g = T_g, T_b = T_b,
        wind = wind, L = L, RH = RH,
        S = exposure_df_hourly[i, "rad_energy"] * (1 - exposure_df_hourly[i, "prop_diffuse_solar"]),
        doy = exposure_df_hourly[i, "doy"], lat = 48.4504,
        lon = -122.9635, hour = exposure_df_hourly[i, "hour"]
      )
    } else {
      Q_rate <- Q_unsteady_noevap( # otherwise, don't allow evaporation
        T_a = T_air, T_g = T_g, T_b = T_b,
        wind = wind, L = L,
        S = exposure_df_hourly[i, "rad_energy"] * (1 - exposure_df_hourly[i, "prop_diffuse_solar"]),
        doy = exposure_df_hourly[i, "doy"], lat = 48.4504,
        lon = -122.9635, hour = exposure_df_hourly[i, "hour"]
      )
    }
    
    #Calculate heat gained or lost and update body temp
    delta_Q <- Q_rate * time_step
    T_b_new <- T_b + delta_Q / C_total
    delta_Tb <- T_b_new - T_b
    delta_Tb <- max(min(delta_Tb, 5), -5)  # cap to ±5°C body temp change in 15 mins
    T_b_new <- T_b + delta_Tb
    
    exposure_df_hourly[i, "body_temp"] <- T_b_new
    exposure_df_hourly[i, "delta_Tb"] <- delta_Tb
    #exposure_df_hourly[i, "Q_stored"] <- NA_real_  # Not tracking this now
    
    T_b <- T_b_new
    
    # Optional: Helpful printing at every timestep to debug
    cat(sprintf(
      "Step %d: evap=%.3e kg/s, vd_body=%.3e, vd_air=%.3e, Δvd=%.3e, fade=%.3f, water_remaining=%.6f kg, Q_rate=%.4f J/s, T_b=%.2f K -> %.2f K, time_step=%ds\n",
      j, evap, vd_body, vd_air, vd_diff, fade_factor, water_remaining, Q_rate, T_b, T_b_new, time_step
    ))
  }
}


# ============================================================
# SECTION 3: Set Tcrit and calculated accumulated degree hours
# ============================================================

#use a distribution of Tcrit rather than a single one
set.seed(100)

thresholds<-rnorm(n=400, mean = 40, sd = 2)

exposure_df_expanded<-expand_grid(
  exposure_df_hourly,
  threshold = thresholds
) %>%
  # Convert body temperature from Kelvin to Celsius
  mutate(body_temp_C = body_temp - 273.15) %>%
  
  # Calculate how much body temp exceeds the threshold; zero if below threshold
  mutate(
    degree_excess = ifelse(body_temp_C > threshold, body_temp_C - threshold, 0),
    degree_hours_step = degree_excess * 0.25
  ) %>%
  
  group_by(length, threshold) %>%
  arrange(exposure_time, .by_group = TRUE) %>%
  mutate(cumulative_degree_hours = cumsum(degree_hours_step)) %>%
  ungroup()

# Summarize the distribution of degree hours
degree_hour_summary <- exposure_df_expanded %>%
  group_by(exposure_time, length) %>%
  summarise(
    dh_mean = mean(cumulative_degree_hours),
    dh_sd = sd(cumulative_degree_hours),
    dh_sem = sd(cumulative_degree_hours) / sqrt(n()),
    dh_q25 = quantile(cumulative_degree_hours, 0.25),
    dh_q75 = quantile(cumulative_degree_hours, 0.75),
    .groups = "drop"
  )

# Rearrange the summary dataframe (helpful for plotting)
degree_hour_summary <- degree_hour_summary %>%
  arrange(length, exposure_time)

# Example heatmap
ggplot(degree_hour_summary, aes(x = length, y = as.numeric(exposure_time) / 3600, fill = dh_mean)) +
  geom_raster(interpolate = TRUE) +
  scale_x_continuous(
    limits = c(0.01, 0.15),
    breaks = seq(0.01, 0.15, by = 0.01)
  ) +
  scale_y_continuous(
    limits = c(2, 9),
    breaks = seq(2, 9, by = 1)
  ) +
  scale_fill_viridis_c(option = "turbo") +
  labs(
    x = "Length (m)",
    y = "Exposure Time (hours)",
    fill = "Mean Degree Hours (°C)"
  ) +

  theme_classic()

