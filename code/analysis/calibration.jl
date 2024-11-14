# ==============================================================
# Main script for calibration of End of Oil Models
# Sets up all input parameters and solves for cost function that equates
# modeled investment to actual investment
# Handles reference case and alternative models / parameters
# Outputs .jld2 files containing all calibrated parameters and
# drilling cost elasticities for all specs
# ==============================================================


# ==============================================================
# Preliminaries --- directories and packages, load model
# ==============================================================
# Get path to root EndOfOil directory
dircurrent = @__DIR__           # Get path to current file
dircurrent = replace(dircurrent, "\\" => "/")
dircurrent = dircurrent * "/"   # Add trailing slash
dirroot = dircurrent[1:end-14]  # Path to root directory

# Paths to output directories
dirfig = dirroot * "output/figures/"
dirfigs = dirroot * "output/figures_slides/"
dirtab = dirroot * "output/tables/"
dirsnt = dirroot * "output/snt/"    # single-number tex files

# Path to intermediate file directory (not in repo)
include(dircurrent * "intfilepath.jl")

# Load packages and models
include(dircurrent * "EOOmodel.jl")
include(dircurrent * "EOOmodel_AKS.jl")



# ==============================================================
# Input parameters -- reference case
# Specifies all but γ, which is solved for below
# ==============================================================
N = 4                       # number of regions
# Initial 2023 production rates (mmbbl/d)
Q0 = 82.757                 # total global crude production
p_q01 = 2.590+9.609+3.251   # core OPEC: Kuwait + SA + UAE
# other OPEC+, incl Azerbaijan, Bahrain, Bruenei, Kazakhstan, Malaysia, 
# Mexico, Oman, Russia, South Sudan, and Sudan. Drop Angola.
p_q02 = 29.985-1.103+0.619+0.201+0.084+1.891+0.508+1.875+1.049+
    10.554+0.148+0.057 - p_q01
p_q04 = 8.42        # US shale oil
p_q03 = Q0 - p_q01 - p_q02 - p_q04  # Rest of world
p_q0 = [p_q01, p_q02, p_q03, p_q04]

# Demand parameters (scalars)
p_delast = -0.5         # demand elasticity at initial production rate and price
p_Pref = 82.49          # initial price, $/bbl
p_dref = Q0             # initial global production rate, mmbbl/d
p_timetozero = 75.0     # years until demand for oil is zero
p_shiftdelin = 0.       # delay in years before demand decline starts
p_drem = 0.             # share of original demand remaining at end of demand decline  

# Initial annual demand growth
p_dgr_ann = 0.0044

# Reserves remaining (billion bbl)
X0 = 1624.           # total global reserves, Rystad 2PCX
p_x01 = 51. + 271. +72.    # core OPEC: Kuwait + SA + UAE
# For other OPEC+, need to infer reserves for Azerbaijan, Bahrain, Brunei, Malaysia,
# Oman, South Sudan, and Sudan based on production
x02pro = 153. * (0.619+0.201+0.084+0.508+1.049+0.148+0.057) /
    (Q0-29.985-12.927-10.554-4.935-4.198-3.402-1.310-1.891-1.875-
    0.281-0.946-1.818-0.391-0.659)
p_x02 = 696-p_x01-13+143+33+23+x02pro
p_x04 = 192 * p_q04 / 12.927  # US shale oil
p_x03 = X0 - p_x01 - p_x02 - p_x04  # Rest of world
p_x0 = [p_x01, p_x02, p_x03, p_x04]

# Annual exponential decline rates in AKS model
p_lambday = [0.08; 0.08; 0.08; 0.3]

# Real annual discount rates
p_r_ann = [0.03, 0.09, 0.09, 0.09]     
# real annual discount rate for emissions
p_Er_ann = p_r_ann[1] 

# Cost intercepts
p_alpha = [5.0, 10.0, 20.0, 30.0]       # MC at q=0, in $/bbl

# Market power
p_mp = true

# Periods per year to simulate
p_T = 2

# Years to simulate
p_maxY = 200

# Tolerance on shadow values
p_mutol = 1e-6

# Combine parameters into reference case touple (excluding γ)
p_ref0 = (q0=p_q0, delast=p_delast, Pref=p_Pref, dref=p_dref, 
    timetozero=p_timetozero, dgr_ann=p_dgr_ann, 
    x0=p_x0, lambday=p_lambday, r_ann=p_r_ann, Er_ann=p_Er_ann, 
    alpha=p_alpha, T=p_T, maxY=p_maxY, mutol=p_mutol, 
    shiftdelin=p_shiftdelin, drem=p_drem, mp = p_mp)



# ==============================================================
# Calibrate γ for non-AKS model
# ==============================================================
gamma_n_guess = [1.6; 2.3; 2.1; 6.0]
m_n_guess = EOO_N(αin=p_alpha, γin=gamma_n_guess, x0in=p_x0, delast=p_delast,
    Pref=p_Pref, dref=p_dref, dgr_ann=p_dgr_ann,
    timetozero=p_timetozero, r_ann=p_r_ann, Er_ann=p_Er_ann, mp=p_mp, T=p_T)
# Solve for γ
p_gamma_n = findgamma(m_n_guess, p_q0, gamma_n_guess)
# Store in touple
p_n = merge(p_ref0, Dict(:gamma => p_gamma_n))
# Get baseline outcomes
m_n = EOO_N(αin=p_n.alpha, γin=p_n.gamma, x0in=p_n.x0, delast=p_n.delast,
    Pref=p_n.Pref, dref=p_n.dref, dgr_ann=p_n.dgr_ann,
    timetozero=p_n.timetozero, r_ann=p_n.r_ann, Er_ann=p_n.Er_ann, mp=p_n.mp, T=p_n.T)
_, pveci_n, qveci_n, mu0i_n, _, _ = optpath(m_n, false, 0., x0(m_n))
# Get extraction cost elasticities
ce_n = mce(m_n,vec(qveci_n[1,:]))
# Get drilling elasticities wrt price
elasttot_n, elast_n = supplyelast(m_n,mu0i_n)


# ==============================================================
# Calibrate γ for reference case model
# ==============================================================
m_guess = EOO_AKS(eoo_n=m_n, lambday=p_lambday, q0=p_q0, pguess=pveci_n)
# Initial guesses of mu0 and gamma
gamma_guess = [1.2, 2.2, 2.0, 5.9]; mu0_guess = [100., 7., 7., 2.]
# Solve for γ
p_gamma = findgamma(m_guess, p_q0, gamma_guess, mu0_guess)
# Store in touple
p_ref = merge(p_n, Dict(:gamma => p_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m = EOO_AKS(eoo_n=m_temp, lambday=p_ref.lambday, q0=p_ref.q0, pguess=pveci_n)
_, _, pveci, dveci, qveci, mu0i, _, _ = 
    optpath(m, q0(m), false, 1, x0(m), mu0_guess)
# Get drilling cost elasticities
ce_ref = mce(m,vec(dveci[1,:]))
# Get drilling elasticities wrt price
elasttot, elast = drillelast(m,q0(m),mu0i,pveci)
ratio_ref = vec(qveci[2,:] ./ qveci[1,:])


# ==============================================================
# Calibrate γ for model in which all resources decline like shale
# ==============================================================
p_dec30_lambday = fill(p_ref.lambday[N], N)
m_guess = EOO_AKS(eoo_n=m_n, lambday=p_dec30_lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.5, 2.3, 2.1, 5.9]; mu0_guess = [24., 2., 2., 2.]
# Solve for γ
p_dec30_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_dec30 = merge(p_ref, Dict(:lambday => p_dec30_lambday, :gamma => p_dec30_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dec30.alpha, γin=p_dec30.gamma, x0in=p_dec30.x0, delast=p_dec30.delast,
    Pref=p_dec30.Pref, dref=p_dec30.dref, dgr_ann=p_dec30.dgr_ann,
    timetozero=p_dec30.timetozero, r_ann=p_dec30.r_ann, Er_ann=p_dec30.Er_ann, mp=p_dec30.mp, T=p_dec30.T)
m_dec30 = EOO_AKS(eoo_n=m_temp, lambday=p_dec30.lambday, q0=p_dec30.q0, pguess=pveci)
_, _, pveci_dec30, dveci_dec30, qveci_dec30, mu0i_dec30, _, _ = 
    optpath(m_dec30, q0(m_dec30), false, 1, x0(m_dec30), mu0_guess)
# Get drilling cost elasticities
ce_dec30 = mce(m_dec30,vec(dveci_dec30[1,:]))
ratio_dec30 = vec(qveci_dec30[2,:] ./ qveci_dec30[1,:])


# ==============================================================
# Calibrate γ for model in which all resources decline like conventional
# ==============================================================
p_dec08_lambday = fill(p_ref.lambday[1], N)
m_guess = EOO_AKS(eoo_n=m_n, lambday=p_dec08_lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.2, 2.2, 2.1, 5.6]; mu0_guess = [100., 7., 7., 7.]
# Solve for γ
p_dec08_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_dec08 = merge(p_ref, Dict(:lambday => p_dec08_lambday, :gamma => p_dec08_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dec08.alpha, γin=p_dec08.gamma, x0in=p_dec08.x0, delast=p_dec08.delast,
    Pref=p_dec08.Pref, dref=p_dec08.dref, dgr_ann=p_dec08.dgr_ann,
    timetozero=p_dec08.timetozero, r_ann=p_dec08.r_ann, Er_ann=p_dec08.Er_ann, mp=p_dec08.mp, T=p_dec08.T)
m_dec08 = EOO_AKS(eoo_n=m_temp, lambday=p_dec08.lambday, q0=p_dec08.q0, pguess=pveci)
_, _, pveci_dec08, dveci_dec08, qveci_dec08, mu0i_dec08, _, _ = 
    optpath(m_dec08, q0(m_dec08), false, 1, x0(m_dec08), mu0_guess)
# Get drilling cost elasticities
ce_dec08 = mce(m_dec08,vec(dveci_dec08[1,:]))
ratio_dec08 = vec(qveci_dec08[2,:] ./ qveci_dec08[1,:])


# ==============================================================
# Calibrate γ for model in which all resources decline at 6% per year
# ==============================================================
p_dec06_lambday = fill(0.06, N)
m_guess = EOO_AKS(eoo_n=m_n, lambday=p_dec06_lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.2, 2.2, 2.1, 5.6]; mu0_guess = [120., 9., 9., 9.]
# Solve for γ
p_dec06_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_dec06 = merge(p_ref, Dict(:lambday => p_dec06_lambday, :gamma => p_dec06_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dec06.alpha, γin=p_dec06.gamma, x0in=p_dec06.x0, delast=p_dec06.delast,
    Pref=p_dec06.Pref, dref=p_dec06.dref, dgr_ann=p_dec06.dgr_ann,
    timetozero=p_dec06.timetozero, r_ann=p_dec06.r_ann, Er_ann=p_dec06.Er_ann, mp=p_dec06.mp, T=p_dec06.T)
m_dec06 = EOO_AKS(eoo_n=m_temp, lambday=p_dec06.lambday, q0=p_dec06.q0, pguess=pveci)
_, _, pveci_dec06, dveci_dec06, qveci_dec06, mu0i_dec06, _, _ = 
    optpath(m_dec06, q0(m_dec06), false, 1, x0(m_dec06), mu0_guess)
# Get drilling cost elasticities
ce_dec06 = mce(m_dec06,vec(dveci_dec06[1,:]))
ratio_dec06 = vec(qveci_dec06[2,:] ./ qveci_dec06[1,:])



# ==============================================================
# Calibrate γ for model with low reserves (Rystad 2PC rather than 2PCX)
# ==============================================================
# First build lowres x0 vector
X0_lowres = 1283.0
p_lowres_x01 = 48+257+69    # core OPEC: Kuwait + SA + UAE
# For other OPEC+, need to infer reserves for Azerbaijan, Bahrain, Brunei, Malaysia,
# Oman, South Sudan, and Sudan based on production
x02pro = 71 * (0.619+0.201+0.084+0.508+1.049+0.148+0.057) /
    (Q0-29.985-12.927-10.554-4.935-4.198-3.402-1.310-1.891-1.875-
    0.281-0.946-1.818-0.391-0.659)
p_lowres_x02 = 638-p_lowres_x01-8+126+27+16+x02pro
p_lowres_x04 = 122 * p_q04 / 12.927  # US shale oil
p_lowres_x03 = X0_lowres - p_lowres_x01 - p_lowres_x02 - p_lowres_x04  # Rest of world
p_lowres_x0 = [p_lowres_x01, p_lowres_x02, p_lowres_x03, p_lowres_x04]
# Instantiate model and calibrate γ
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_lowres_x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [0.7, 2.1, 1.8, 5.3]; mu0_guess = [140., 15., 15., 10.]
# Solve for γ
p_lowres_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_lowres = merge(p_ref, Dict(:x0 => p_lowres_x0, :gamma => p_lowres_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_lowres.alpha, γin=p_lowres.gamma, x0in=p_lowres.x0, delast=p_lowres.delast,
    Pref=p_lowres.Pref, dref=p_lowres.dref, dgr_ann=p_lowres.dgr_ann,
    timetozero=p_lowres.timetozero, r_ann=p_lowres.r_ann, Er_ann=p_lowres.Er_ann, mp=p_lowres.mp, T=p_lowres.T)
m_lowres = EOO_AKS(eoo_n=m_temp, lambday=p_lowres.lambday, q0=p_lowres.q0, pguess=pveci)
_, _, pveci_lowres, dveci_lowres, qveci_lowres, mu0i_lowres, _, _ = 
    optpath(m_lowres, q0(m_lowres), false, 1, x0(m_lowres), mu0_guess)
# Get drilling cost elasticities
ce_lowres = mce(m_lowres,vec(dveci_lowres[1,:]))
ratio_lowres = vec(qveci_lowres[2,:] ./ qveci_lowres[1,:])


# ==============================================================
# Calibrate γ for model with high reserves (IEA rather than Rystad 2PCX)
# ==============================================================
X0_highres = 2602.0
p_highres_x0 = p_x0 * X0_highres / X0
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_highres_x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [2.3, 2.3, 2.1, 6.1]; mu0_guess = [30., 1., 1., 1.]
# Solve for γ
p_highres_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_highres = merge(p_ref, Dict(:x0 => p_highres_x0, :gamma => p_highres_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_highres.alpha, γin=p_highres.gamma, x0in=p_highres.x0, delast=p_highres.delast,
    Pref=p_highres.Pref, dref=p_highres.dref, dgr_ann=p_highres.dgr_ann,
    timetozero=p_highres.timetozero, r_ann=p_highres.r_ann, Er_ann=p_highres.Er_ann, mp=p_highres.mp, T=p_highres.T)
m_highres = EOO_AKS(eoo_n=m_temp, lambday=p_highres.lambday, q0=p_highres.q0, pguess=pveci)
_, _, pveci_highres, dveci_highres, qveci_highres, mu0i_highres, _, _ = 
    optpath(m_highres, q0(m_highres), false, 1, x0(m_highres), mu0_guess)
# Get drilling cost elasticities
ce_highres = mce(m_highres,vec(dveci_highres[1,:]))
ratio_highres = vec(qveci_highres[2,:] ./ qveci_highres[1,:])


# ==============================================================
# Calibrate γ for models in which all types discount at alternative rates
# ==============================================================
# Automate with function
function findgammadisc(p_ref, d, gamma_guess, mu0_guess, N, pveci)
    p_d_r_ann = fill(d, N)
    m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
        Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
        timetozero=p_ref.timetozero, r_ann=p_d_r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
    m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
    # Solve for γ
    p_d_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
    # Store in touple
    p_d = merge(p_ref, Dict(:r_ann => p_d_r_ann, :gamma => p_d_gamma))
    # Get baseline outcomes
    m_temp = EOO_N(αin=p_d.alpha, γin=p_d.gamma, x0in=p_d.x0, delast=p_d.delast,
        Pref=p_d.Pref, dref=p_d.dref, dgr_ann=p_d.dgr_ann,
        timetozero=p_d.timetozero, r_ann=p_d.r_ann, Er_ann=p_d.Er_ann, mp=p_d.mp, T=p_d.T)
    m_d = EOO_AKS(eoo_n=m_temp, lambday=p_d.lambday, q0=p_d.q0, pguess=pveci)
    _, _, pveci_d, dveci_d, qveci_d, mu0i_d, _, _ = 
        optpath(m_d, q0(m_d), false, 1, x0(m_d), mu0_guess)
    # Get drilling cost elasticities
    ce_d = mce(m_d,vec(dveci_d[1,:]))
    ratio_d = vec(qveci_d[2,:] ./ qveci_d[1,:])
    return p_d, ce_d, ratio_d, mu0i_d, pveci_d
end
# Start with 3% discounting and work up from there
# Initial guesses of mu0 and gamma
gamma_guess = [1.4, 1.3, 1.2, 3.0]; mu0_guess = [100., 120., 120., 30.]
p_dall03, ce_dall03, ratio_dall03, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.03, gamma_guess, mu0_guess, N, pveci)
p_dall04, ce_dall04, ratio_dall04, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.04, p_dall03.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall05, ce_dall05, ratio_dall05, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.05, p_dall04.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall06, ce_dall06, ratio_dall06, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.06, p_dall05.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall07, ce_dall07, ratio_dall07, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.07, p_dall06.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall08, ce_dall08, ratio_dall08, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.08, p_dall07.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall09, ce_dall09, ratio_dall09, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.09, p_dall08.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall10, ce_dall10, ratio_dall10, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.10, p_dall09.gamma.+0.1, mu0i_d/1.5, N, pveci_d)
p_dall12, ce_dall12, ratio_dall12, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.12, p_dall10.gamma.+0.1, mu0i_d/1.3, N, pveci_d)
p_dall14, ce_dall14, ratio_dall14, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.14, p_dall12.gamma.+0.1, mu0i_d/1.3, N, pveci_d)
p_dall16, ce_dall16, ratio_dall16, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.16, p_dall14.gamma.+0.1, mu0i_d/1.3, N, pveci_d)
p_dall18, ce_dall18, ratio_dall18, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.18, p_dall16.gamma.+0.05, mu0i_d/1.1, N, pveci_d)
p_dall20, ce_dall20, ratio_dall20, mu0i_d, pveci_d = 
    findgammadisc(p_ref, 0.20, p_dall18.gamma.+0.05, mu0i_d/1.1, N, pveci_d)


    
# ==============================================================
# Calibrate γ for non-AKS model, with low reserves
# ==============================================================
p_n_lowres_x0 = p_lowres_x0
# Initial guess of gamma
gamma_n_lowres_guess = [0.7, 2.1, 1.8, 5.3]
m_n_lowres_guess = EOO_N(αin=p_ref.alpha, γin=gamma_n_guess, x0in=p_n_lowres_x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
# Solve for γ
p_gamma_n_lowres = findgamma(m_n_lowres_guess, p_ref.q0, gamma_n_lowres_guess)
# Store in touple
p_n_lowres = merge(p_ref0, Dict(:x0 => p_n_lowres_x0, :gamma => p_gamma_n_lowres))
# Get baseline outcomes
m_n_lowres = EOO_N(αin=p_n_lowres.alpha, γin=p_n_lowres.gamma, x0in=p_n_lowres.x0, delast=p_n_lowres.delast,
    Pref=p_n_lowres.Pref, dref=p_n_lowres.dref, dgr_ann=p_n_lowres.dgr_ann,
    timetozero=p_n_lowres.timetozero, r_ann=p_n_lowres.r_ann, Er_ann=p_n_lowres.Er_ann, mp=p_ref.mp, T=p_n_lowres.T)
_, pveci_n_lowres, qveci_n_lowres, mu0i_n_lowres, _, _ = optpath(m_n_lowres, false, 0., x0(m_n_lowres))
# Get extraction cost elasticities
ce_n_lowres = mce(m_n_lowres,vec(qveci_n_lowres[1,:]))
# Get drilling elasticities wrt price
elasttot_n_lowres, elast_n_lowres = supplyelast(m_n_lowres,mu0i_n_lowres)



# ==============================================================
# Calibrate γ for non-AKS model, with 3% discounting
# ==============================================================
p_n_dall03_r_ann = fill(0.03, N)
# Initial guess of gamma
gamma_n_dall03_guess = [1.4, 1.3, 1.2, 3.0]
m_n_dall03_guess = EOO_N(αin=p_ref.alpha, γin=gamma_n_guess, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_n_dall03_r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
# Solve for γ
p_gamma_n_dall03 = findgamma(m_n_dall03_guess, p_ref.q0, gamma_n_dall03_guess)
# Store in touple
p_n_dall03 = merge(p_ref0, Dict(:r_ann => p_n_dall03_r_ann, :gamma => p_gamma_n_dall03))
# Get baseline outcomes
m_n_dall03 = EOO_N(αin=p_n_dall03.alpha, γin=p_n_dall03.gamma, x0in=p_n_dall03.x0, delast=p_n_dall03.delast,
    Pref=p_n_dall03.Pref, dref=p_n_dall03.dref, dgr_ann=p_n_dall03.dgr_ann,
    timetozero=p_n_dall03.timetozero, r_ann=p_n_dall03.r_ann, Er_ann=p_n_dall03.Er_ann, mp=p_ref.mp, T=p_n_dall03.T)
_, pveci_n_dall03, qveci_n_dall03, mu0i_n_dall03, _, _ = optpath(m_n_dall03, false, 0., x0(m_n_dall03))
# Get extraction cost elasticities
ce_n_dall03 = mce(m_n_dall03,vec(qveci_n_dall03[1,:]))
# Get drilling elasticities wrt price
elasttot_n_dall03, elast_n_dall03 = supplyelast(m_n_dall03,mu0i_n_dall03)



# ==============================================================
# Set up models with alternative demand declines
# ==============================================================
p_T50 = merge(p_ref, Dict(:timetozero => 50))
p_T100 = merge(p_ref, Dict(:timetozero => 100))
p_del05 = merge(p_ref, Dict(:shiftdelin => 05))
p_del05_dall09 = merge(p_dall09, Dict(:shiftdelin => 05))
p_dem85 = merge(p_ref, Dict(:drem => 0.15, :maxY => 300, :mutol => 1e-10))



# ==============================================================
# Calibrate γ for model with larger core OPEC and 9% discounting
# ==============================================================
# Move Russia to core OPEC
q0_R = 10.554
x0_R = 143
p_R_q0 = [p_ref.q0[1]+q0_R, p_ref.q0[2]-q0_R, p_ref.q0[3], p_ref.q0[4]] 
p_R_x0 = [p_ref.x0[1]+x0_R, p_ref.x0[2]-x0_R, p_ref.x0[3], p_ref.x0[4]]
p_R_r_ann = fill(0.09, N)
# Instantiate model and calibrate γ
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_R_x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_R_r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_R_q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [0.87, 3.4, 2.0, 5.9] ; mu0_guess = [5., 5., 5., 2.]
# Solve for γ
p_R_gamma = findgamma(m_guess, p_R_q0, gamma_guess, mu0_guess)
# Store in touple
p_R = merge(p_ref, Dict(:q0 => p_R_q0, :x0 => p_R_x0, :r_ann => p_R_r_ann, :gamma => p_R_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_R.alpha, γin=p_R.gamma, x0in=p_R.x0, delast=p_R.delast,
    Pref=p_R.Pref, dref=p_R.dref, dgr_ann=p_R.dgr_ann,
    timetozero=p_R.timetozero, r_ann=p_R.r_ann, Er_ann=p_R.Er_ann, mp=p_R.mp, T=p_R.T)
m_R = EOO_AKS(eoo_n=m_temp, lambday=p_R.lambday, q0=p_R.q0, pguess=pveci)
_, _, pveci_R, dveci_R, qveci_R, mu0i_R, _, _ = 
    optpath(m_R, q0(m_R), false, 1, x0(m_R), mu0_guess)
# Get drilling cost elasticities
ce_R = mce(m_R,vec(dveci_R[1,:]))
ratio_R = vec(qveci_R[2,:] ./ qveci_R[1,:])



# ==============================================================
# Calibrate γ for model with no market power
# ==============================================================
p_nmp_mp = false
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_nmp_mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.5, 2.2, 2.0, 5.9]; mu0_guess = [110., 7., 7., 2.]
# Solve for γ
p_nmp_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_nmp = merge(p_ref, Dict(:mp => p_nmp_mp, :gamma => p_nmp_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_nmp.alpha, γin=p_nmp.gamma, x0in=p_nmp.x0, delast=p_nmp.delast,
    Pref=p_nmp.Pref, dref=p_nmp.dref, dgr_ann=p_nmp.dgr_ann,
    timetozero=p_nmp.timetozero, r_ann=p_nmp.r_ann, Er_ann=p_nmp.Er_ann, mp=p_nmp.mp, T=p_nmp.T)
m_nmp = EOO_AKS(eoo_n=m_temp, lambday=p_nmp.lambday, q0=p_nmp.q0, pguess=pveci)
_, _, pveci_nmp, dveci_nmp, qveci_nmp, mu0i_nmp, _, _ = 
    optpath(m_nmp, q0(m_nmp), false, 1, x0(m_nmp), mu0_guess)
# Get drilling cost elasticities
ce_nmp = mce(m_nmp,vec(dveci_nmp[1,:]))
ratio_nmp = vec(qveci_nmp[2,:] ./ qveci_nmp[1,:])


# ==============================================================
# Calibrate γ for model with high demand elasticity
# ==============================================================
p_delhigh_delast = -0.6
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_delhigh_delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.7, 2.2, 2.0, 5.9]; mu0_guess = [90., 7., 7., 2.]
# Solve for γ
p_delhigh_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_delhigh = merge(p_ref, Dict(:delast => p_delhigh_delast, :gamma => p_delhigh_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_delhigh.alpha, γin=p_delhigh.gamma, x0in=p_delhigh.x0, delast=p_delhigh.delast,
    Pref=p_delhigh.Pref, dref=p_delhigh.dref, dgr_ann=p_delhigh.dgr_ann,
    timetozero=p_delhigh.timetozero, r_ann=p_delhigh.r_ann, Er_ann=p_delhigh.Er_ann, mp=p_delhigh.mp, T=p_delhigh.T)
m_delhigh = EOO_AKS(eoo_n=m_temp, lambday=p_delhigh.lambday, q0=p_delhigh.q0, pguess=pveci)
_, _, pveci_delhigh, dveci_delhigh, qveci_delhigh, mu0i_delhigh, _, _ = 
    optpath(m_delhigh, q0(m_delhigh), false, 1, x0(m_delhigh), mu0_guess)
# Get drilling cost elasticities
ce_delhigh = mce(m_delhigh,vec(dveci_delhigh[1,:]))
ratio_delhigh = vec(qveci_delhigh[2,:] ./ qveci_delhigh[1,:])



# ==============================================================
# Calibrate γ for model with low demand elasticity
# Do this manually since the γ search causes non-convergence
# ==============================================================
p_dellow_delast = -0.4; p_dellow_gamma = [0.548, 2.2017, 2.0485, 5.892]
mu0_guess = [120., 9., 7., 2.5]
# Store in touple
p_dellow = merge(p_ref, Dict(:delast => p_dellow_delast, :gamma => p_dellow_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dellow.alpha, γin=p_dellow.gamma, x0in=p_dellow.x0, delast=p_dellow.delast,
    Pref=p_dellow.Pref, dref=p_dellow.dref, dgr_ann=p_dellow.dgr_ann,
    timetozero=p_dellow.timetozero, r_ann=p_dellow.r_ann, Er_ann=p_dellow.Er_ann, mp=p_dellow.mp, T=p_dellow.T)
m_dellow_guess = EOO_AKS(eoo_n=m_temp, lambday=p_dellow.lambday, q0=p_dellow.q0, pguess=pveci)
# Get better initial price vector
_, _, pveci_dellow_guess, _, _, _, _, _ = 
    optpath(m_dellow_guess, q0(m_dellow_guess), false, 1, x0(m_dellow_guess), mu0_guess)
# Re-do model with new intiial price vec
m_dellow = EOO_AKS(eoo_n=m_temp, lambday=p_dellow.lambday, q0=p_dellow.q0, pguess=pveci_dellow_guess)
_, _, pveci_dellow, dveci_dellow, qveci_dellow, mu0i_dellow, _, _ = 
    optpath(m_dellow, q0(m_dellow), false, 1, x0(m_dellow), mu0_guess)
# Get drilling cost elasticities
ce_dellow = mce(m_dellow,vec(dveci_dellow[1,:]))
ratio_dellow = vec(qveci_dellow[2,:] ./ qveci_dellow[1,:])



# ==============================================================
# Calibrate γ for model with no demand growth
# ==============================================================
p_dng_dgr_ann = 0.0
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_dng_dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.4, 2.3, 2.1, 6.0]; mu0_guess = [95., 6., 6., 2.]
# Solve for γ
p_dng_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_dng = merge(p_ref, Dict(:dgr_ann => p_dng_dgr_ann, :gamma => p_dng_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dng.alpha, γin=p_dng.gamma, x0in=p_dng.x0, delast=p_dng.delast,
    Pref=p_dng.Pref, dref=p_dng.dref, dgr_ann=p_dng.dgr_ann,
    timetozero=p_dng.timetozero, r_ann=p_dng.r_ann, Er_ann=p_dng.Er_ann, mp=p_dng.mp, T=p_dng.T)
m_dng = EOO_AKS(eoo_n=m_temp, lambday=p_dng.lambday, q0=p_dng.q0, pguess=pveci)
_, _, pveci_dng, dveci_dng, qveci_dng, mu0i_dng, _, _ = 
    optpath(m_dng, q0(m_dng), false, 1, x0(m_dng), mu0_guess)
# Get drilling cost elasticities
ce_dng = mce(m_dng,vec(dveci_dng[1,:]))
ratio_dng = vec(qveci_dng[2,:] ./ qveci_dng[1,:])



# ==============================================================
# Calibrate γ for model with high (OPEC forecast) demand growth
# ==============================================================
p_dhg_dgr_ann = 0.015
m_n_guess = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_dhg_dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [1.0, 2.1, 1.9, 5.8]; mu0_guess = [105., 8., 8., 3.]
# Solve for γ
p_dhg_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_dhg = merge(p_ref, Dict(:dgr_ann => p_dhg_dgr_ann, :gamma => p_dhg_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_dhg.alpha, γin=p_dhg.gamma, x0in=p_dhg.x0, delast=p_dhg.delast,
    Pref=p_dhg.Pref, dref=p_dhg.dref, dgr_ann=p_dhg.dgr_ann,
    timetozero=p_dhg.timetozero, r_ann=p_dhg.r_ann, Er_ann=p_dhg.Er_ann, mp=p_dhg.mp, T=p_dhg.T)
m_dhg = EOO_AKS(eoo_n=m_temp, lambday=p_dhg.lambday, q0=p_dhg.q0, pguess=pveci)
_, _, pveci_dhg, dveci_dhg, qveci_dhg, mu0i_dhg, _, _ = 
    optpath(m_dhg, q0(m_dhg), false, 1, x0(m_dhg), mu0_guess)
# Get drilling cost elasticities
ce_dhg = mce(m_dhg,vec(dveci_dhg[1,:]))
ratio_dhg = vec(qveci_dhg[2,:] ./ qveci_dhg[1,:])



# ==============================================================
# Calibrate γ for model with higher cost function intercepts
# ==============================================================
p_higha_alpha = [15.0, 20.0, 30.0, 40.0]       # MC at q=0, in $/bbl
m_n_guess = EOO_N(αin=p_higha_alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [0.7, 1.9, 1.7, 4.8]; mu0_guess = [100., 7., 7., 2.]
# Solve for γ
p_higha_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_higha = merge(p_ref, Dict(:alpha => p_higha_alpha, :gamma => p_higha_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_higha.alpha, γin=p_higha.gamma, x0in=p_higha.x0, delast=p_higha.delast,
    Pref=p_higha.Pref, dref=p_higha.dref, dgr_ann=p_higha.dgr_ann,
    timetozero=p_higha.timetozero, r_ann=p_higha.r_ann, Er_ann=p_higha.Er_ann, mp=p_higha.mp, T=p_higha.T)
m_higha = EOO_AKS(eoo_n=m_temp, lambday=p_higha.lambday, q0=p_higha.q0, pguess=pveci)
_, _, pveci_higha, dveci_higha, qveci_higha, mu0i_higha, _, _ = 
    optpath(m_higha, q0(m_higha), false, 1, x0(m_higha), mu0_guess)
# Get drilling cost elasticities
ce_higha = mce(m_higha,vec(dveci_higha[1,:]))
ratio_higha = vec(qveci_higha[2,:] ./ qveci_higha[1,:])



# ==============================================================
# Calibrate γ for model with even higher cost function intercepts
# ==============================================================
p_highera_alpha = [20.0, 25.0, 35.0, 45.0]       # MC at q=0, in $/bbl
m_n_guess = EOO_N(αin=p_highera_alpha, γin=p_ref.gamma, x0in=p_ref.x0, delast=p_ref.delast,
    Pref=p_ref.Pref, dref=p_ref.dref, dgr_ann=p_ref.dgr_ann,
    timetozero=p_ref.timetozero, r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=p_ref.mp, T=p_ref.T)
m_guess = EOO_AKS(eoo_n=m_n_guess, lambday = p_ref.lambday, q0=p_ref.q0, pguess=pveci)
# Initial guesses of mu0 and gamma
gamma_guess = [0.462, 1.745, 1.555, 4.1935]; mu0_guess = [100., 7., 7., 2.]
# Solve for γ
p_highera_gamma = findgamma(m_guess, p_ref.q0, gamma_guess, mu0_guess)
# Store in touple
p_highera = merge(p_ref, Dict(:alpha => p_highera_alpha, :gamma => p_highera_gamma))
# Get baseline outcomes
m_temp = EOO_N(αin=p_highera.alpha, γin=p_highera.gamma, x0in=p_highera.x0, delast=p_highera.delast,
    Pref=p_highera.Pref, dref=p_highera.dref, dgr_ann=p_highera.dgr_ann,
    timetozero=p_highera.timetozero, r_ann=p_highera.r_ann, Er_ann=p_highera.Er_ann, mp=p_highera.mp, T=p_highera.T)
m_highera = EOO_AKS(eoo_n=m_temp, lambday=p_highera.lambday, q0=p_highera.q0, pguess=pveci)
_, _, pveci_highera, dveci_highera, qveci_highera, mu0i_highera, _, _ = 
    optpath(m_highera, q0(m_highera), false, 1, x0(m_highera), mu0_guess)
# Get drilling cost elasticities
ce_highera = mce(m_highera,vec(dveci_highera[1,:]))
ratio_highera = vec(qveci_highera[2,:] ./ qveci_highera[1,:])



# ==============================================================
# Export touples and elasticities to file
# ==============================================================
jldsave(dirint * "cal.jld2"; p_ref0, p_n, p_ref, p_dec30, p_dec08, p_dec06,
    p_lowres, p_highres, p_dall03, p_dall04, p_dall05, p_dall06, p_dall07, p_dall08,
    p_dall09, p_dall10, p_dall12, p_dall14, p_dall16, p_dall18, p_dall20,
    p_n_lowres, p_n_dall03, p_T50, p_T100, p_del05, p_del05_dall09, p_dem85, p_R, p_nmp,
    p_delhigh, p_dellow, p_dng, p_dhg, p_higha, p_highera)
jldsave(dirint * "elast.jld2"; elasttot, elast, ce_n, ce_ref, ce_dec30, ce_dec08, ce_dec06,
    ce_lowres, ce_highres, ce_dall03, ce_dall04, ce_dall05, ce_dall06, ce_dall07, ce_dall08,
    ce_dall09, ce_dall10, ce_dall12, ce_dall14, ce_dall16, ce_dall18, ce_dall20,
    ce_n_lowres, ce_n_dall03, ce_R, ce_nmp, ce_delhigh, ce_dellow, ce_dng, ce_dhg, ce_higha, ce_highera)

