# ==============================================================
# Script for generating results given custom inputs
# Accompanies "The End of Oil" by Ryan Kellogg
# ==============================================================


# ==============================================================
# PRELIMINARIES --- set directories and load model functions
# No need to for users to modify this section
# ==============================================================
# Get path to root EndOfOil directory
dircurrent = @__DIR__           # Get path to current file
dircurrent = replace(dircurrent, "\\" => "/")
dircurrent = dircurrent * "/"   # Add trailing slash

# Load packages and models
include(dircurrent * "EOOmodel.jl")
include(dircurrent * "EOOmodel_AKS.jl")
include(dircurrent * "PlotFunctions.jl")



# ==============================================================
# INPUT PARAMETERS
# All are currently set to match the reference case as described in the paper
# The model can accommodate fewer or more than the four regions in the paper. 
# If the number of regions is changed, then all input vectors must be updated 
# to have the same number of elements.
# ==============================================================

# Flag for whether to include AKS (2018) investment dynamics in the model
# Default is true. Change to false to model extraction without investment
p_AKS = true

# Flag for whether to re-calibrate the cost function slope parameters gamma.
# Default is false. Change to true to re-calibrate the cost function slopes
# "True" is recommended if you change other parameters. If changes are made and
# this is left as "false", the model may not match 2023 actual investment
# and production in the baseline scenario.
p_estslopes = false

# Flag for whether to run the unanticipated decline scenario. Default is true.
# This scenario substantially increases computation time, since the model must
# be run for each period of the simulation. Change to false to skip this scenario.
p_unant = true

# Vector of initial 2023 production rates (mmbbl/d)
p_q0 = [15.45, 30.418, 28.469, 8.42]

# Demand parameters (scalars)
p_delast = -0.5         # demand elasticity at initial production rate and price
p_Pref = 82.49          # initial price, $/bbl
p_dref = sum(p_q0)      # initial global production rate, mmbbl/d
p_timetozero = 75.0     # years until demand for oil is zero in decline scenario
p_shiftdelin = 0.       # delay in years before demand decline starts
p_drem = 0.             # share of original demand remaining at end of demand decline  
p_dgr_ann = 0.0044      # initial annual demand growth

# Vector of reserves remaining (billion bbl)
p_x0 = [394.0, 541.7769281476599, 563.1638933886595, 125.05917846368067] 

# Annual exponential decline rates in AKS model
# Ignored if p_AKS = false
p_lambday = [0.08; 0.08; 0.08; 0.3]

# Real annual discount rates
p_r_ann = [0.03, 0.09, 0.09, 0.09]     
# real annual discount rate for computing present discounted emissions
p_Er_ann = 0.03

# Drilling MC function intercepts, in $/bbl
p_alpha = [5.0, 10.0, 20.0, 30.0]       # MC at q=0, in $/bbl

# Drilling MC function slopes.
# Units are $/bbl per mmbbl/d invested for models with AKS investment dynamics
# Units are $/bbl per mmbbl/d of production for models without investment dynamics
# These are used as calibration starting values if p_estslopes = true. See README for 
# information on good guesses for these values depending on other input parameters.
p_gamma = [30.3947, 53.8906, 50.0906, 36.2644]

# Flag for whether core OPEC (resource 1) exerts market power. "true" means yes.
p_mp = true

# Initial guesses for initial shadow values in baseline scenario, in $/bbl
# See README for information on good guesses for these values depending on 
# other input parameters
mu0_bbl_guess = [30.8797, 3.59828, 2.7732, 2.37266]

# Simulation length / tolerance parameters. Users likely do not need to change these
p_T = 2             # Periods per year to simulate
p_maxY = 200        # Years to simulate
p_mutol = 1e-6      # Tolerance on shadow values



# ==============================================================
# RUN MODEL. NO NEED FOR USERS TO EDIT THIS CODE
# ==============================================================
# If AKS model, convert units of gamma from $/bbl to mmbbl/d of investment per period
# to $/bbl per mmbbl/d of steady-state production
if p_AKS
    lambdap = 1 .- (1 .- p_lambday).^(1/p_T)     # per-period decline rate
    p_gammain = p_gamma .* lambdap
else
    p_gammain = p_gamma
end
# Build tuple of all input parameters (not including the three flags)
p_all = (q0=p_q0, delast=p_delast, Pref=p_Pref, dref=p_dref, 
    timetozero=p_timetozero, dgr_ann=p_dgr_ann, 
    x0=p_x0, lambday=p_lambday, r_ann=p_r_ann, Er_ann=p_Er_ann, 
    alpha=p_alpha, gamma=p_gammain, T=p_T, maxY=p_maxY, mutol=p_mutol, 
    shiftdelin=p_shiftdelin, drem=p_drem, mp = p_mp)
# Run model
Baseline, Ant, Unant, PctInc = runall(p_all, p_AKS, p_estslopes, p_unant, mu0_bbl_guess)
# Plot combined drilling, production, and prices
if p_unant==true; L = 3; else; L = 2; end   # set number of scenarios to plot
combined_plot = comboplot(Baseline, Ant, Unant, p_T, true, L)
# Region plots of production. Store as vector of plots (each element is a region)
region_plots = []
for i in 1:length(p_q0)
    region_plot_i = regionplot(Baseline.qvec, Ant.qvec, Unant.qvec, p_T, i, true, true, true, L)
    push!(region_plots, region_plot_i)
end



# ==============================================================
# DICTIONARY OF MODEL OUTPUTS
# ==============================================================
# The outputs Baseline, Ant, and Unant are tuples that have the same structure of fields.
# Baseline corresponds to the baseline demand scenario, Ant to the anticipated demand 
#   decline scenario, and Unant to the unanticipated decline scenario.
# Let N denote the number of resources and TY the number of periods simulated
# The fields of these tuples are described below. They can be accessed by typing, e.g., `Baseline.Q`:
    # D: Nx1 vector of cumulative drilling investment, in mmbbl/d of capacity added for each resource
    # Q: Nx1 vector of cumulative extraction, in billion bbl for each resource
    # pvec: TYx1 vector of oil prices in each period, in $/bbl
    # dvec: NxTY matrix of per-period drilling investment for each resource, in mmbbl/d per period
    # qvec: NxTY matrix of per-period extraction for each resource, in mmbbl/d
    # muvec: NxTY matrix of current shadow values for each resource, in $ per mmbbl/d
    # mu0: Nx1 vector of initial shadow values, in $ per mmbbl/d
    # mu0_bbl: Nx1 vector of initial shadow values, in $/bbl
    # Te: Nx1 vector of period in which there is strictly positive drilling, for each resource
    # mc_bbl: Nx1 vector of initial marginal investment cost for each resource, in $/bbl
    # QT: scalar of total cumulative extraction across all resources, in billion bbl
    # PVQT: scalar of present value of total cumulative extraction across all resources, in billion bbl
    # PVQ: Nx1 vector of present value of cumulative extraction for each resource, in billion bbl
    # pivec: NxTY matrix of per-period profits for each resource, in $billion
    # PVpi: Nx1 vector of present value of profits for each resource, in $billion per period
# Note that if p_unant==false, the fields of Unant will be placeholder values of 0 or 1.

# The output PctInc contains the following scalar fields that summarize percentage changes in 
#   cumulative outcomes, aggregated across all resource types:
    # PctIncDdDiTot: cumulative investment in anticipated decline relative to baseline
    # PctIncDdDmTot: cumulative investment in anticipated decline relative to unanticipated decline
    # PctIncQdQiTot: cumulative extraction in anticipated decline relative to baseline
    # PctIncQdQmTot: cumulative extraction in anticipated decline relative to unanticipated decline
    # PctIncQmQiTot: cumulative extraction in unanticipated decline relative to baseline
    # PctIncPVQdQiTot: present discounted cumulative extraction in anticipated decline relative to baseline
    # PctIncPVQdQmTot: present discounted cumulative extraction in anticipated decline relative to unanticipated decline
    # PctIncPVQmQiTot: present discounted cumulative extraction in unanticipated decline relative to baseline

# Plots
# combined_plot is a vertically stacked three-panel plot of drilling investment, production, and prices,
#   aggregated across all resource types
# region_plots is a vector of plots, each of which shows production under each scenario (baseline, anticipated decline,
#   and unanticipated decline) for a single resource type. 
#   The first element of region_plots corresponds to the first resource type, and so on.
#   The elements of region_plots can be accessed by typing region_plots[i], where i is the index of the resource type.
