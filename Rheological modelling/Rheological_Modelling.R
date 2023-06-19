## Script to calculate deformation mechanisms in quartz, calcite, chlorite, and illite
## Used for rheological modelling of Asanzaki data
## Created by Madison Frank, updated 2023/06/07

library(plotly)
library(SciViews)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
## General parameters

# molar gas constant (J K-1 mol-1)
R <- 8.31446261815324 

# temperature from Raman analysis
Temp_C <- 286 # in C
Temp_K <- Temp_C + 273.15 # in K

# Calculate depth (km) in subduction zone based on geothermal gradient of 15 C/km (Penniston-Dorland et al., 2015)
depth <- Temp_C / 15

# gravity (km/s)
gravity <- 9.8 * 10^-3 

# Solid density of qtz (kg/m3) from Robertson 1988
qtz_density <- 2648 

# Overburden (MPa) calculated from qtz density 
overburden <- qtz_density * depth * gravity

# Fluid density as a function of depth and Temp_C (kg/m3) calculated @ 494.79 MPa w/NIST Lemmon et al.
fluid_density <- 1002.3

# Provide differential stress to calculate strain rate
diff_stress <- 10^seq(-2, 5, by = 0.1) # MPa
shear_stress <- diff_stress/sqrt(3) # convert to shear stress (Paterson & Olgaard, 2000))

# generic bulk strain rates of 'gouge' to model over - odd stepsize is to ensure same length as frictional shear (161)
strain_rates <- 10^seq(-28, 1, by = 0.41) 

# provide generic velocities (s^-1) to model over - odd stepsize to match length of diff_stress (71)
slip_vel <- 10^seq(-28, 1, by = 0.41)

# Pore fluid factors
pore_fluid_factor <- c(0.4, 0.6, 0.8, 0.9, 0.95, 0.999)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Calculate true outcrop thickness from average dip and width measured in field

# MDM
MDM_dip <- 46.3 # dip (corrected to radians)
MDM_width <- 35 + 8 # measured from Google Earth across outcrop between Bay 2/3 parallel to dip direction + cliff height 
MDM_thickness <- sin(MDM_dip * pi/180)*MDM_width # true thickness of MDM

# BLM
BLM_dip <- 44 # dip (corrected to radians)
BLM_width <- 26.1 # measured farthest SW BLM outcrop
BLM_thickness <- sin(BLM_dip * pi/180)*BLM_width # true thickness of BLM


# ILLITE #--------------------------------------------------------------------------------------------------------------------------------------------------------

## Illite parameters

# frictional coefficient (Hartog et al., 2012)
il_u_0 <- 0.4 

# rate-dependence of friction with time (Hartog et al., 2012)
il_du_dlnV_T <- -0.005 

# initial velocity converted from 1 um/s
il_V_0 <- 10^-6 

# rate dependent friction coefficient (Condit 2019, eq. 3) - NOTE log is actually ln in R
il_u_VT <- il_u_0 + (il_du_dlnV_T * log(slip_vel/il_V_0))

#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Coulomb Frictional slip (Condit et al., 2019, Hartog et al., 2012)

# create an empty dataframe to include frictional slip at each pore fluid factor
il_stress_DF <- data.frame(Velocity = slip_vel, u = il_u_VT)

# loop through each pore fluid factor
for (i in 1:length(pore_fluid_factor)) {

  # calculate shear stress
  t_il <- il_u_VT * ((1 - pore_fluid_factor[i]) * overburden)

  # add results to dataframe
  il_stress_DF <- cbind(il_stress_DF, t_il)
}

# rename dataframe columns
colnames(il_stress_DF) <- c("Velocity", "u", paste("Pf_factor", pore_fluid_factor, sep = "_"))

#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot illite stress results

il_stress_fig <- plot_ly(il_stress_DF, x = ~Velocity) %>%

  add_trace(y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "Pf = 0.4", line = list(color = "darkblue")) %>%
  add_trace(y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "Pf = 0.6", line = list(color = "blue")) %>%
  add_trace(y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "Pf = 0.8", line = list(color = "dodgerblue")) %>%
  add_trace(y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "Pf = 0.9", line = list(color = "deepskyblue")) %>%
  add_trace(y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "Pf = 0.95", line = list(color = "cyan")) %>%
  add_trace(y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "Pf = 0.999", line = list(color = "paleturquoise")) %>%

  layout(title = "Illite Rheology",
         xaxis = list(title = "Slip Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log"),
         yaxis = list(title = "Shear Stress (MPa)", type = "log", showexponent = "all", exponentformat = "power"),
         font = list(size = 13),
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = FALSE
  )

il_stress_fig


# QUARTZ #--------------------------------------------------------------------------------------------------------------------------------------------------------

## Quartz parameters

# Average quartz grain size (um)
d_qtz <- 40 * 10^-6

# Molar volume of quartz from Berman 1988 (m3/mol)
V_m <- 2.269*10^-5

# phenomenological rate coefficient for dissolution (Rimstidt & Barnes, 1980)
k_sp <- 10^(1.174 + (-2.028*10^-3 * Temp_K) + (-4158 / Temp_K)) 

#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Quartz dislocation creep (Hirth et al., 2001)

# Material constant (MPa^-n / s)
qtz_A_dcl <- 10^-11.2 

# activation energy (J/mol)
qtz_Q_dcl <- 135*1000 

# fugacity at overburden and Temp_C (MPa) from T. Withers calculator
f_h2o <- 0.0463195437568053  * 10^3

# Dislocation strain rate
qtz_disloc_strain <- qtz_A_dcl * f_h2o^1 * diff_stress^4 * exp(-qtz_Q_dcl / (R*Temp_K))


#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Diffusion controlled dissolution-precipitation solution - Bos/Spiers 2002a model (eq. 17)

# aspect ratio of grains
B <- 2.2

# angle of microlithon (rads) - The leading edge of the grain is inclined at angle a to the horizontal
alpha <- 23 * pi/180 

# factor expressing proportion of grain boundary area undergoing sliding
P <- 0.75 


# friction component (results of illite sliding above x proportion of active foliae)
qtz_fric_stress <- P * il_stress_DF[-c(1,2)]


# Phenomenological rate coefficient for dissolution - eq from Rimstidt & Barnes 1980
K_sp <- (R * Temp_K * d_qtz) / (B^2 * k_sp * V_m)

# Pressure solution component
qtz_PS_stress <- (overburden * tan(alpha) * K_sp * strain_rates) / ((overburden * tan(alpha)) + (K_sp * strain_rates))


# create empty dataframe to add combined friction + pressure solution for all pore fluid factors
qtz_fric_PS_shear_stress_DF <- data.frame(Velocity = slip_vel)

# loop through pore fluid factors
for (i in 1:length(pore_fluid_factor)) {
  
  # adjust overburden for current pore fluid pressure
  norm_stress <- ((1 - pore_fluid_factor[i]) * overburden)
  
  # calculate pressure solution component for current pore fluid pressure
  temp_PS_stress <- (norm_stress * tan(alpha) * K_sp * strain_rates) / ((norm_stress * tan(alpha)) + (K_sp * strain_rates))
  
  # add frictional sliding and pressure solution components 
  qtz_fric_PS_shear_stress <- qtz_fric_stress[i] + temp_PS_stress
    
  # add to dataframe for storing combined fric + pressure solution
  qtz_fric_PS_shear_stress_DF <- cbind(qtz_fric_PS_shear_stress_DF, qtz_fric_PS_shear_stress)
}

# Rename columns to match pore fluid factor
colnames(qtz_fric_PS_shear_stress_DF) <- c("Velocity", paste("Pf_factor", pore_fluid_factor, sep = "_"))


#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot quartz stress-strain results

# Create master dataframe with all qtz deformation mechanisms vs shear stress - exclude microphysical model as it will be plotted 
# again velocity (x) rather than shear stress (y)
qtz_strain_DF <- data.frame(ShearStress = diff_stress/sqrt(3), # convert differential stress to shear stress (Paterson & Olgaard)
                            Qtz_DislocationCreep = qtz_disloc_strain
                            )



qtz_strain_fig <- plot_ly(qtz_strain_DF, y = ~ShearStress) %>%  
  
  add_trace(x = ~Qtz_DislocationCreep, mode = "line", type = "scatter", name = "Dislocation Creep") %>%

  add_trace(data = qtz_fric_PS_shear_stress_DF, x = ~Velocity, y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "BS Microphysical 0.4 (Bos & Spiers 2002)") %>%
  add_trace(data = qtz_fric_PS_shear_stress_DF, x = ~Velocity, y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "BS Microphysical 0.8 (Bos & Spiers 2002)") %>%
  add_trace(data = qtz_fric_PS_shear_stress_DF, x = ~Velocity, y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "BS Microphysical 0.9 (Bos & Spiers 2002)") %>%

  layout(title = "Quartz Rheology",
         
         xaxis = list(title = "Strain Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log", 
                      range = list(log10(10^-18), log10(10^-4))),
         
         yaxis = list(title = "Shear Stress (MPa)", type = "log", showexponent = "all", exponentformat = "power", 
                      range = list(log10(10^-2), log10(10^4))),
         
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h')
         )

qtz_strain_fig



# CHLORITE #--------------------------------------------------------------------------------------------------------------------------------------------------------

## Chlorite parameters

# frictional coefficient (Okamoto Fig 7)
chl_u_0 <- 0.27 

# rate-dependence (Okamoto Fig 8)
chl_du_dlnV_T <- -0.001 

# (converted from 1 um/s)
chl_V_0 <- 10^-6 

# rate dependent friction coefficient (Condit 2019, eq. 3) - NOTE log is actually ln in R
chl_u_VT <- chl_u_0 + (chl_du_dlnV_T * log(slip_vel/chl_V_0)) 


#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Coulomb Frictional slip (Condit et al., 2019, Okamoto et al., 2019)

# create with dataframe to include frictional slip at each pore fluid factor
chl_stress_DF <- data.frame(Velocity = slip_vel, u = chl_u_VT)

# loop through each pore fluid factor
for (i in 1:length(pore_fluid_factor)) {
  
  # calculate shear stress 
  t_chlorite <- chl_u_VT * ((1 - pore_fluid_factor[i]) * overburden)
  
  # add results to dataframe
  chl_stress_DF <- cbind(chl_stress_DF, t_chlorite)
}

# rename columns
colnames(chl_stress_DF) <- c("Velocity", "u", paste("Pf_factor", pore_fluid_factor, sep = "_"))


#--------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot chlorite stress results

chl_stress_fig <- plot_ly(chl_stress_DF, x = ~Velocity) %>%  
  
  add_trace(y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "Pf = 0.4", line = list(color = "darkblue")) %>%
  add_trace(y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "Pf = 0.6", line = list(color = "blue")) %>%
  add_trace(y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "Pf = 0.8", line = list(color = "dodgerblue")) %>%
  add_trace(y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "Pf = 0.9", line = list(color = "deepskyblue")) %>%
  add_trace(y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "Pf = 0.95", line = list(color = "cyan")) %>%
  add_trace(y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "Pf = 0.999", line = list(color = "paleturquoise")) %>%
  
  layout(title = "Chlorite Rheology",
         xaxis = list(title = "Slip Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log"),
         yaxis = list(title = "Shear Stress (MPa)", type = "log", showexponent = "all", exponentformat = "power"),
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h')  )

chl_stress_fig





# CALCITE #--------------------------------------------------------------------------------------------------------------------------------------------------------
## Calcite parameters

# grain size diameter (um) of microcrystalline calcite calculated from MTEX - NOTE Walker 1990 specifically states grain size in um

d_Ca <- 7.7707


#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Calcite dislocation creep > 1 s (Rutter 1974)

# activation energy (J/mol) - Reading et al. (1984), Renner et al. (2002)
ca_Q_dcl <- 200*1000 

# frequency/pre-exponential material constant (s^-1) - French & Condit, 2019
ca_A_dcl <- 10^6 

# concentration (MPa^-1) - French & Condit, 2019
ca_C <- 0.035 

# calculate dislocation creep
Ca_disloc_strain <- ca_A_dcl * exp(-ca_Q_dcl / (R*Temp_K)) * exp(ca_C * diff_stress)

# remove values at rates < 1 s
Ca_disloc_strain[diff_stress < 1] <- NA

# Cut off asymptotic projection for plotting
Ca_disloc_strain[Ca_disloc_strain > 1] <- NA


#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Calcite grain boundary sliding diffusion creep (Walker 1990)  

# activation energy (J/mol)
ca_Q_diff <- 190*1000 

# frequency/pre-exponential material constant (s^-1)
ca_A_diff_low <- 10^4.9 # differential stress < 25 MPa
ca_A_diff_high <- 100 # differential stress > 25 MPa

# calculate diffusion strain rate at differential stress < 25 MPa
Ca_diff_strain_low <- ca_A_diff_low * exp(-(ca_Q_diff) / (R*Temp_K)) * diff_stress^1.7 * d_Ca^-1.9 

# remove values at differential stress > 25
Ca_diff_strain_low[diff_stress > 25] <- NA

# calculate diffusion strain rate at differential stress > 25 MPa
Ca_diff_strain_high <- ca_A_diff_high * exp(-(ca_Q_diff) / (R*Temp_K)) * diff_stress^3.3 * d_Ca^-1.3 

# remove values at differential stress < 25
Ca_diff_strain_high[diff_stress < 25] <- NA



#--------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot calcite stress results

# create master dataframe with all calcite deformation mechanisms
Ca_strain_DF <- data.frame(ShearStress = diff_stress/sqrt(3),
                           Ca_DislocationCreep = Ca_disloc_strain,
                           Ca_DiffusionCreep_Low = Ca_diff_strain_low,
                           Ca_DiffusionCreep_High = Ca_diff_strain_high)


Ca_strain_fig <- plot_ly(Ca_strain_DF, y = ~ShearStress) %>%  
  
  add_trace(x = ~Ca_DislocationCreep, mode = "lines", type = "scatter", name = "Dislocation Creep (>80 MPa)", line = list(color = '#2ca02c')) %>%
  add_trace(x = ~Ca_DiffusionCreep_Low, mode = "lines", type = "scatter", name = "Diffusion Creep (<25 MPa)", line = list(color = '#ff7f0e')) %>%  
  add_trace(x = ~Ca_DiffusionCreep_High, mode = "lines", type = "scatter", name = "Diffusion Creep (>25 MPa)", line = list(color = '#d62728')) %>%
  
  layout(title = "Calcite Rheology",
         xaxis = list(title = "Strain Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log"),
         yaxis = list(title = "Shear Stress (MPa)", type = "log", showexponent = "all", exponentformat = "power"),
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h')
         )


Ca_strain_fig



### STRAIN RATE MODEL ###---------------------------------------------------------------------------------------------------------------------------------------------------
## MDM

# rename dataframes to avoid risk of overwriting
MDM_strain_DF_Pf <- qtz_strain_DF#[c(1,2)]
MDM_stress_DF_Pf <- qtz_fric_PS_shear_stress_DF

# Plot all MDM deformation mechanisms as stress vs strain
MDM_strain_Pf_fig <- plot_ly(MDM_strain_DF_Pf, y = ~ShearStress) %>%  
  
  add_trace(x = ~Qtz_DislocationCreep, mode = "lines", type = "scatter", name = "Dislocation Creep", line = list(color = "#2ca02c")) %>%

  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "BS Pf = 0.4", line = list(color = "darkorange")) %>%
  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "BS Pf = 0.6", line = list(color = "orange")) %>%
  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "BS Pf = 0.8", line = list(color = "darkgoldenrod")) %>%
  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "BS Pf = 0.9", line = list(color = "gold")) %>%
  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "BS Pf = 0.95", line = list(color = "khaki")) %>%
  add_trace(data = MDM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "BS Pf = 0.999", line = list(color = "lemonchiffon")) %>%
  
  layout(title = paste0("MDM Strain Rate"),
         
         xaxis = list(title = "Strain Rate (s<sup>-1</sup>)", exponentformat = "power", 
                      type = "log", range = list(log10(10^-24), log10(10^-4))),
         
         yaxis = list(title = "Shear Stress (MPa)", showexponent = "all", exponentformat = "power",
                      type = "log", range = list(log10(10^-2), log10(10^4))),
         
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h')
         )

MDM_strain_Pf_fig





#--------------------------------------------------------------------------------------------------------------------------------------------------------
## BLM


# rename dataframes to avoid risk of overwriting
BLM_strain_DF_Pf <- Ca_strain_DF
BLM_stress_DF_Pf <- chl_stress_DF

# plot all BLM deformation mechanisms as stress vs strain
BLM_strain_fig <- plot_ly(BLM_strain_DF_Pf, y = ~ShearStress) %>%  

  add_trace(x = ~Ca_DislocationCreep, mode = "line", type = "scatter", name = "Dislocation Creep (>80 MPa)") %>%
  add_trace(x = ~Ca_DiffusionCreep_Low, mode = "line", type = "scatter", name = "Diffusion Creep (<25 MPa)") %>%  
  add_trace(x = ~Ca_DiffusionCreep_High, mode = "line", type = "scatter", name = "Diffusion Creep (>25 MPa)") %>%
  
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "Pf = 0.4", line = list(color = "darkblue")) %>%
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "Pf = 0.6", line = list(color = "blue")) %>%
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "Pf = 0.8", line = list(color = "dodgerblue")) %>%
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "Pf = 0.9", line = list(color = "deepskyblue")) %>%
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "Pf = 0.95", line = list(color = "cyan")) %>%
  add_trace(data = BLM_stress_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "Pf = 0.999", line = list(color = "paleturquoise")) %>%

  layout(title = paste0("BLM Strain Rate"),
         
         xaxis = list(title = "Strain Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log", 
                      range = list(log10(10^-16), log10(10^-4))),
         
         yaxis = list(title = "Shear Stress (MPa)", showexponent = "all", exponentformat = "power", type = "log", 
                      range = list(log10(10^-2), log10(10^4))),
         
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h')
         )


BLM_strain_fig

### SLIP RATE MODEL ###--------------------------------------------------------------------------------------------------------------------------------------------
## MDM

# multiple thickness by 10 to have shear zone thickness between 100 - 350m (Rowe et al., 2013)
MDM_assumed_thickness <- MDM_thickness * 10

# rename dataframes to avoid overwrite
MDM_slipRate_DF_Pf <- MDM_strain_DF_Pf

# Convert strain rate to slip rate
MDM_slipRate_DF_Pf[2] <- MDM_slipRate_DF_Pf[2] * MDM_assumed_thickness 

# rename dataframes to avoid overwrite
MDM_fric_slipRate_DF_Pf <- MDM_stress_DF_Pf

# convert velocity to slip rate
MDM_fric_slipRate_DF_Pf[1] <- MDM_fric_slipRate_DF_Pf[1] * MDM_assumed_thickness


# plot all MDM deformation mechanisms as shear stress vs slip rate
MDM_slipRate_fig <- plot_ly(MDM_slipRate_DF_Pf, y = ~ShearStress) %>%  
  
  add_trace(x = ~Qtz_DislocationCreep, mode = "lines", type = "scatter", name = "Dislocation Creep", line = list(color = "#2ca02c")) %>%

  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "Pf = 0.4", line = list(color = "darkorange")) %>%
  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "Pf = 0.6", line = list(color = "orange")) %>%
  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "Pf = 0.8", line = list(color = "darkgoldenrod")) %>%
  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "Pf = 0.9", line = list(color = "gold")) %>%
  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "Pf = 0.95", line = list(color = "khaki")) %>%
  add_trace(data = MDM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "Pf = 0.999", line = list(color = "lemonchiffon")) %>%
  
  layout(title = paste0("MDM Slip Rate @ ", round(MDM_thickness, 0), "m thickness"),
         
         xaxis = list(title = "Slip Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log", 
                      range = list(log10(10^-22), log10(10^-4))),
         
         yaxis = list(title = "Shear Stress (MPa)", showexponent = "all", exponentformat = "power", type = "log", 
                      range = list(log10(10^-2), log10(10^4))),
         
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h'),
         
         # add rectangles who show aseismic and slow slip rates
         shapes = list(
           list(type = "rect",
                fillcolor = "blue", line = list(color = "blue"), opacity = 0.2, xref = "x",
                x0 = 3.2*10^-10, x1 = 3.2*10^-9,#10^-11, x1 = 10^-10
                y0 = MDM_strain_DF_Pf[1,1], y1 = MDM_strain_DF_Pf[nrow(MDM_strain_DF_Pf),1], yref = "y"),
           
           list(type = "rect",
                fillcolor = "green", line = list(color = "green"), opacity = 0.2, xref = "x",
                x0 = 1.15*10^-8, x1 = 1.15*10^-7,
                y0 = MDM_strain_DF_Pf[1,1], y1 = MDM_strain_DF_Pf[nrow(MDM_strain_DF_Pf),1], yref = "y"))
  )

MDM_slipRate_fig





#--------------------------------------------------------------------------------------------------------------------------------------------------------
## BLM

# multiple thickness by 10 to have shear zone thickness between 100 - 350m (Rowe et al., 2013)
BLM_assumed_thickness <- BLM_thickness * 10

# create new dataframe where all strain rate columns are converted to slip rate, leave stress (col 1) as is
BLM_slipRate_DF_Pf <- data.frame(BLM_strain_DF_Pf[1], BLM_strain_DF_Pf[-1] * BLM_assumed_thickness) # Convert to slip rate 

# rename dataframes to avoid overwrite
BLM_fric_slipRate_DF_Pf <- BLM_stress_DF_Pf

# Convert velocity to slip rate
BLM_fric_slipRate_DF_Pf[1] <- BLM_fric_slipRate_DF_Pf[1] * BLM_assumed_thickness


# plot all MDM deformation mechanisms as shear stress vs slip rate
BLM_slipRate_fig <- plot_ly(BLM_slipRate_DF_Pf, y = ~ShearStress) %>%  
  
  add_trace(x = ~Ca_DislocationCreep, mode = "line", type = "scatter", name = "Dislocation Creep (>80 MPa)") %>%

  add_trace(x = ~Ca_DiffusionCreep_Low, mode = "line", type = "scatter", name = "Diffusion Creep (<25 MPa)") %>%
  add_trace(x = ~Ca_DiffusionCreep_High, mode = "line", type = "scatter", name = "Diffusion Creep (>25 MPa)") %>%

  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.4, mode = "lines", type = "scatter", name = "Pf = 0.4", line = list(color = "darkblue")) %>%
  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.6, mode = "lines", type = "scatter", name = "Pf = 0.6", line = list(color = "blue")) %>%
  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.8, mode = "lines", type = "scatter", name = "Pf = 0.8", line = list(color = "dodgerblue")) %>%
  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.9, mode = "lines", type = "scatter", name = "Pf = 0.9", line = list(color = "deepskyblue")) %>%
  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.95, mode = "lines", type = "scatter", name = "Pf = 0.95", line = list(color = "cyan")) %>%
  add_trace(data = BLM_fric_slipRate_DF_Pf, x = ~Velocity, y = ~Pf_factor_0.999, mode = "lines", type = "scatter", name = "Pf = 0.999", line = list(color = "paleturquoise")) %>%
  
  layout(title = paste0("BLM Slip Rate @ ", round(BLM_thickness, 0), "m thickness"),
         
         xaxis = list(title = "Slip Rate (s<sup>-1</sup>)", exponentformat = "power", type = "log", 
                      range = list(log10(10^-14), log10(10^-4))),
         
         yaxis = list(title = "Shear Stress (MPa)", showexponent = "all", exponentformat = "power", type = "log", 
                      range = list(log10(10^-2), log10(10^4))),
         
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10),
         showlegend = TRUE,
         legend = list(orientation = 'h'),
         
         # add rectangles who show aseismic and slow slip rates
         shapes = list(
           list(type = "rect",
                fillcolor = "blue", line = list(color = "blue"), opacity = 0.2, xref = "x",
                x0 = 3.2*10^-10, x1 = 3.2*10^-9,#10^-11, x1 = 10^-10
                y0 = BLM_strain_DF_Pf[1,1], y1 = BLM_strain_DF_Pf[nrow(BLM_strain_DF_Pf),1], yref = "y"),
           list(type = "rect",
                fillcolor = "green", line = list(color = "green"), opacity = 0.2, xref = "x",
                x0 = 1.15*10^-8, x1 = 1.15*10^-7,
                y0 = BLM_strain_DF_Pf[1,1], y1 = BLM_strain_DF_Pf[nrow(BLM_strain_DF_Pf),1], yref = "y"),
           
           # add rectangle to show shear stress from Ca twins 
           list(type = "rect",
                fillcolor = "grey", line = list(color = "grey"), opacity = 0.2, xref = "x",
                x0 = BLM_strain_DF_Pf[1,3], x1 = BLM_strain_DF_Pf[nrow(BLM_strain_DF_Pf),1],
                y0 = 97/2, y1 = 252/2, yref = "y"))
  )


BLM_slipRate_fig

