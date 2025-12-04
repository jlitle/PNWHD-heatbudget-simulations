# =============================================================================
# Script name: heatbudget_functions.R
# Author: Jack Litle
# Date: 2025-12-02
# Purpose: Define functions necessary for running biophysical simulations 
# ============================================================================

#Load required packages
library(TrenchR)

#Define functions for biophysical simulations
initial_water<- function(L){
  water<-0.37*191*(L^3.53)
  ifelse(water<=0.0001563, 0.0001563, water) #minimum starting water
}

shell_mass<-function(L){
  shellmass<-57.6*L^2.92
  ifelse(shellmass<=0.001208, 0.001208, shellmass) #minimum shell mass
}

tissue_mass<-function(L){
  tissuemass<-191*L^3.53
  ifelse(tissuemass<=0.0004224367, 0.0004224367, tissuemass) #minimum tissue mass
}

h_m<-function(wind,L){
  (2*26e-6)/L + (0.60*(((wind*0.66*L)/(16e-6))^0.5)*(16/26.5)^(1/3)*26e-6)/L
}

area<- function(L){
  raw_area<-0.1*((1.08*(L^2) + 0.0461*L - 0.0016))  #area open to evaporation
  ifelse(raw_area < 0.1*0.0002275, 0.0002275, raw_area)
}

area_mussel<- function(L){ 
  raw_area<- ((1.08*(L^2) + 0.0461*L - 0.0016)) #mussel area as a function of length
  ifelse(raw_area < 0.0002275, 0.0002275, raw_area)
}


SatVapDense<-function(T){ 
  T=T-273.15
  (5.018 +0.32321*T +(8.1847e-3*T^2)+(3.1243e-4*T^3))/1000 #saturation vapor density at temperature T in Celsius
}

evap_rate<- function(wind,L, vd_body, vd_air){
  (h_m(wind=wind, L=L)*area(L)*(vd_body-vd_air)) 
}


Re<- function(wind, L){ #Reynold's number
  (wind*0.66*L)/16e-6  #kinematic viscosity of air in m2/s
}

hc<- function(wind, L){ #heat transfer coefficient
  (0.67*(Re(wind=wind, L=L)^0.42)*0.026)/0.66*L #0.026 is the thermal conductivity of air in W/M2/K
}

solar_elev_angle<-function(doy, lat, lon, hour){
  (90 -(zenith_angle(doy=doy, lat=lat, lon = lon, hour = hour)))*(pi/180)
}

projected_area<-function(theta, a, b){ #theta is solar elevation angle in radians
  pi*a*b*sqrt((cos(theta)^2 + (b/a)^2 * (sin(theta)^2)))
}

# ====================================
# Unsteady Heat Transfer (Q) Functions
# ====================================

#old function
sol_unsteady<-function(doy, lat, lon, hour, S, L){
  0.75*(sin(solar_elev_angle(doy=doy, lat=lat, lon=lon, hour=hour)))^-1 * 0.15*area_mussel(L=L)* S
}

rad_sky_unsteady<- function(T_a, T_b, L) {
  4*5.67E-8*0.85^0.75*0.5*area_mussel(L=L)*T_a^3*(T_b-0.85^0.25*T_a)
}

rad_ground_unsteady<- function(T_g, T_b, L) {
  4*5.67E-8*0.5*area_mussel(L=L)*T_g^3*(T_b-T_g)
}

cond_unsteady<- function(T_g, T_b, L) {
  0.6*(0.5*0.5*L)^-1*0.05*area_mussel(L=L)*(T_b-T_g)
}

conv_unsteady<- function(T_a,T_b,wind,L) {
  hc(wind=wind, L=L)*area_mussel(L=L)*(T_b-T_a)
}

evap_unsteady<- function(wind, L, T_b, T_a, RH) {
  2.48E6*evap_rate(wind=wind, L=L, vd_body = SatVapDense(T_b), vd_air = RH * (SatVapDense(T_a))) #2.48 = latent heat of vaporization of water in J/kg
}

#full unsteady heat budgets
Q_unsteady_evap<-function(T_a, T_g, T_b, wind, L, RH, S, doy, lat, lon, hour){ #with evaporation
  sol_unsteady(doy = doy, lat = lat, lon = lon, hour = hour, S = S, L = L) - 
    rad_sky_unsteady(T_a = T_a, T_b = T_b, L=L) - 
    rad_ground_unsteady(T_g = T_g, T_b = T_b, L = L) - 
    cond_unsteady(T_g = T_g, T_b = T_b, L = L) -
    conv_unsteady(T_a = T_a, T_b = T_b, wind = wind, L=L) - 
    evap_unsteady(wind = wind, L=L, T_b = T_b, T_a = T_a, RH = RH )
}

Q_unsteady_noevap<-function(T_a, T_g, T_b, wind, L, S, doy, lat, lon, hour){ #without evaporation
  sol_unsteady(doy = doy, lat = lat, lon = lon, hour = hour, S = S, L = L) - 
    rad_sky_unsteady(T_a = T_a, T_b = T_b, L=L) - 
    rad_ground_unsteady(T_g = T_g, T_b = T_b, L = L) - 
    cond_unsteady(T_g = T_g, T_b = T_b, L = L) -
    conv_unsteady(T_a = T_a, T_b = T_b, wind = wind, L=L)
}
