# =============================================================================
# Script name: summarize_model_outputs.R
# Author: Jack Litle
# Date: 2025-12-02
# Purpose: Summarize results from model outputs for PNWHD mussel heatwave study
# Notes: Only core model summarization pipeline included; figure/statistical scripts not included
# =============================================================================

#Load required packages
library(ggplot2)
library(TrenchR)
library(tidyverse)
library(splines)


# =========================================================================
# SECTION 1: Convert exposure times to tidal heights using lookup summary files
# and extract summary information from model outputs
# =========================================================================

jun28lookup<-read.csv("new_june28_lookup.csv", sep = ",", header = TRUE) # each day has its own lookup day to convert exposure time to tidal elevation

jun28output<-read.csv("ctmaxdist_june28_dhsummary.csv", sep = ",", header = TRUE) # loading daily model outputs generated using run_transfer_simulations.R

summary_sheet<-rbind(jun25lookup, jun26lookup, jun27lookup, jun28lookup) # combine all daily lookup files from PNWHD

output_df<-rbind(jun25output, jun26output, jun27output, jun28output) # combine all daily output files from PNWHD

# Now extract summary info from large combined output_df

# Ensure output_df is sorted
output_df <- output_df[order(output_df$doy, output_df$length, output_df$exposure_time), ]

# Preallocate result column
summary_sheet$degree_hours <- NA_real_
summary_sheet$dh_sd <- NA_real_
summary_sheet$dh_sem <- NA_real_
summary_sheet$dh_q25 <- NA_real_
summary_sheet$dh_q75 <- NA_real_

# Loop over each unique day
for (d in unique(summary_sheet$doy)) {
  
  # Then for each mussel length on that day
  for (m in unique(summary_sheet$length)) {
    
    # Get all rows in summary_sheet that match this day and length
    summary_rows <- which(summary_sheet$doy == d & summary_sheet$length == m)
    
    # Get all matching rows from output_df for this group
    group_df <- output_df[output_df$doy == d & output_df$length == m, ]
    
    # Skip if group is empty
    if (nrow(group_df) == 0) next
    
    # Loop over all summary rows in this group
    for (i in summary_rows) {
      which_row <- summary_sheet$which.row[i]
      
      # Check if the row exists in group_df
      if (nrow(group_df) >= which_row) {
        summary_sheet$degree_hours[i] <- group_df$dh_mean[which_row]
        summary_sheet$exposure_time_from_other[i]<-(group_df$exposure_time[which_row])/60
        summary_sheet$dh_sd[i] <- group_df$dh_sd[which_row]
        summary_sheet$dh_sem[i] <- group_df$dh_sem[which_row]
        summary_sheet$dh_q25[i] <- group_df$dh_q25[which_row]
        summary_sheet$dh_q75[i] <- group_df$dh_q75[which_row]
      }
    }
  }
}

write.csv(summary_sheet, "ctmaxdist_tide_height_full_output.csv", row.names = FALSE)

# =========================================================================
# SECTION 2: Track mean accumulated degree hours (and variance)
# =========================================================================

# propogating standard error of the mean
error_prop_fun<-function(x){ #SE adds quadratically
  sqrt(sum(x^2))
}

# summarizing accumulated degree hours and sem for plotting
heatmap_df <- summary_sheet %>%
  group_by(length, tidal_height) %>%
  summarize(
    degree_hours = sum(degree_hours, na.rm = TRUE),
    dh_sem = error_prop_fun(dh_sem),
    .groups = "drop"
  )

view(summary_sheet) 

#df conducive to plotting
write.csv(heatmap_df, "heatmap_df_ctmaxdist_tide_height_full_output.csv", row.names = FALSE)


# =========================================================================
# SECTION 3: Modelling desiccation probability as a function of tidal elevation
# and size
# =========================================================================

# Load water loss summary
water_loss<-read.csv("water_loss_full_output.csv", sep = ",", header = TRUE)

# Organize df by tidal elevation and compute number desiccated
waterloss_df <- water_loss %>%
  group_by(tidal_height) %>%
  summarise(
    num_dry = sum(water_lost >= 25),
    num_total = n(),
    probs = sum(water_lost >= 25) / n(),
    .groups = "drop"
  )


#modelling desiccation prob as a function of tidal height grouping sizes
#need to model tidal elevation using splines due to nonlinearity

library(splines)
model <- glm(cbind(num_dry, num_total - num_dry) ~ ns(tidal_height, df = 3),
             family = binomial(link = "logit"),
             data = waterloss_df)

summary(model)

#modelling the tidal height x size interaction

model_indiv <- glm(
  I(water_lost >= 25) ~ ns(tidal_height, df = 3) * length,
  family = binomial(link = "logit"),
  data = water_loss
)
summary(model_indiv)

#anova to just summarize the main effects

library(car)
Anova(model_indiv, type = "II") 

# To visualize the fit with confidence intervals:

summary_df <- waterloss_df %>%
  mutate(
    fitted_prob = predict(model, type = "response"),
    se = predict(model, type = "link", se.fit = TRUE)$se.fit,
    fit_link = predict(model, type = "link"),
    lower_link = fit_link - 1.96 * se,
    upper_link = fit_link + 1.96 * se,
    lower_prob = plogis(lower_link),
    upper_prob = plogis(upper_link)
  )

#example plot
ggplot(summary_df, aes(x = tidal_height)) +
  geom_point(aes(y = num_dry / num_total), size = 3.5, color = "black") +
  geom_line(aes(y = fitted_prob), size = 3, color = "black") +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.4, fill = "grey") +
  scale_y_continuous(
    limits = c(0.3, 0.71),
    breaks = seq(0.3, 0.7, by = 0.1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(-0.21, 1.31),
    breaks = seq(-0.2, 1.4, by = 0.1),
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 1.9),   
    axis.ticks = element_line(size = 1.9),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

#Modelling using 3 representative sizes as in Figure 3 of main text
water_loss<-read.csv("water_loss_full_output.csv", sep = ",", header = TRUE)


# Create binary drying column
water_loss$dry <- as.integer(water_loss$water_lost >= 25)

#summary of raw probabilities if you want it 
obs_summary <- water_loss %>%
  group_by(tidal_height, length) %>%
  summarise(prob_dry = mean(dry), n = n())

# Fit GLM with tidal height spline + size
# Scale size
water_loss$length_s <- scale(water_loss$length)
water_loss$length_s <- as.numeric(scale(water_loss$length))

# Additive model with reduced spline flexibility
model_add <- glm(
  dry ~ ns(tidal_height, df = 2) + length_s,
  family = binomial(link = "logit"),
  data = water_loss
)

#Alternatively, use an interaction model
model_int <- glm(
  dry ~ ns(tidal_height, df = 2) * length_s,
  family = binomial(link = "logit"),
  data = water_loss
)

# Test interaction
anova(model_add, model_int, test = "Chisq") #likelihood ratio says interaction matters

#nov4fig
sizes_raw <- c(0.05, 0.08, 0.10)
sizes_scaled <- (sizes_raw - mean(water_loss$length)) / sd(water_loss$length)
sizes_scaled
sizes_scaled <- as.numeric((sizes_raw - mean(water_loss$length)) / sd(water_loss$length))
newdat <- expand.grid(
  tidal_height = seq(min(water_loss$tidal_height),
                     max(water_loss$tidal_height),
                     length.out = 200),
  length_s = sizes_scaled
)

p <- predict(model_add, newdata = newdat, type = "link", se.fit = TRUE)

newdat$fit_link <- p$fit
newdat$se_link  <- p$se.fit

# 95% CI on link scale
newdat$fit_lo <- newdat$fit_link - 1.96 * newdat$se_link
newdat$fit_hi <- newdat$fit_link + 1.96 * newdat$se_link

# convert link → probability
newdat$pred   <- plogis(newdat$fit_link)
newdat$pred_lo <- plogis(newdat$fit_lo)
newdat$pred_hi <- plogis(newdat$fit_hi)

newdat$size_class <- factor(
  newdat$length_s,
  levels = sizes_scaled,
  labels = c("Small (0.05)", "Medium (0.08)", "Large (0.10)")
)

#example plot
ggplot(newdat, aes(x = tidal_height, y = pred, color = size_class, fill = size_class)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.20, color = NA) +
  geom_line(size = 1.3) +
  scale_color_manual(values = c(
    "#1b9e77",  # green
    "#d95f02",  # orange
    "#7570b3",  # purple
    '#7570b3'
    # Add or remove colors as needed
  )) +
  labs(x = "Tidal height",
       y = "Probability of drying",
       color = "Mussel size",
       fill = "Mussel size") +
  theme_classic(base_size = 14)

ggplot(newdat, aes(x = tidal_height, y = pred, color = size_class)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = size_class), alpha = 0.2, color = NA) +
  geom_line(size = 3) +
  scale_color_manual(values = c(
    "#1b9e77",  # green
    "#d95f02",  # orange
    "#7570b3"   # purple
  )) +
  scale_fill_manual(values = c(
    "#1b9e77",  # green
    "#d95f02",  # orange
    "#7570b3"   # purple
  )) +
  labs(x = "Tidal height",
       y = "Probability of drying",
       color = "Mussel size",
       fill = "Mussel size") +
  theme_classic(base_size = 14)

ggplot(newdat, aes(x = tidal_height, y = pred, color = size_class)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = size_class), alpha = 0.2, color = NA) +
  geom_line(size = 3) +
  scale_color_manual(values = c(
    "#1b9e77",  # green
    "#d95f02",  # orange
    "#7570b3"   # purple
  )) +
  scale_fill_manual(values = c(
    "#1b9e77",  # green
    "#d95f02",  # orange
    "#7570b3"   # purple
  )) +
  labs(x = NULL, y = NULL, color = NULL, fill = NULL) +
  theme_classic(base_size = 14) +

  theme(legend.position = "none")
