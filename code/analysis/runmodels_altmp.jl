# ==============================================================
# Script for running set of End of Oil Models
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
# Load input parameters
# ==============================================================
data = load(dirint * "cal.jld2")


# ==============================================================
# Run alternative reserves
# ==============================================================
results_altmp = Dict()
Threads.@threads for spec in ["p_ref", "p_R", "p_nmp"]
    p = data[spec]
    results_altmp[spec] = runAKSmodel(p)
end
Oi_ref, Od_ref, Om_ref, Op_ref = results_altmp["p_ref"]
Oi_R, Od_R, Om_R, Op_R = results_altmp["p_R"]
Oi_nmp, Od_nmp, Om_nmp, Op_nmp = results_altmp["p_nmp"]


# ==============================================================
# Market power goes away during decline --- custom run
# ==============================================================
Oi_mpg = Oi_ref; Om_mpg = Om_ref    # Baseline and unanticipated decline same as ref case
# Define new model with market power going away
p_ref = data["p_ref"]
p_mpg = merge(p_ref, Dict(:mp => false))  # No market power
m_mpg_temp = EOO_N(αin=p_mpg.alpha, γin=p_mpg.gamma, x0in=p_mpg.x0, 
    delast=p_mpg.delast, Pref=p_mpg.Pref, dref=p_mpg.dref, 
    dgr_ann=p_mpg.dgr_ann, timetozero=p_mpg.timetozero, 
    r_ann=p_mpg.r_ann, Er_ann=p_mpg.Er_ann, mp=p_mpg.mp, T=p_mpg.T,
    shiftdelin=p_mpg.shiftdelin, drem=p_mpg.drem, maxY=p_mpg.maxY, mutol=p_mpg.mutol)
m_mpg = EOO_AKS(eoo_n=m_mpg_temp, lambday=p_mpg.lambday, q0=p_mpg.q0, pguess=Oi_ref.pvec)
# Anticipated demand decline
Dd, Qd, pvecd, dvecd, qvecd, mu0d, muvecd, Ted = optpath(m_mpg, q0(m_mpg), true, 1, x0(m_mpg))
QTd, PVQTd, PVQd, PVpid, pivecd = sum_Q_pi(m_mpg, qvecd, pvecd, dvecd)
Od_mpg = (D = Dd, Q = Qd, pvec = pvecd, dvec = dvecd, qvec = qvecd,
        mu0 = mu0d, muvec = muvecd, Te = Ted, QT = QTd, 
        PVQT = PVQTd, PVQ = PVQd, PVpi = PVpid, pivec = pivecd)
# Percent changes
PctIncDdDiTot = sum(Dd-Oi_mpg.D) / sum(Oi_mpg.D) * 100; PctIncDdDmTot = sum(Dd-Om_mpg.D) / sum(Om_mpg.D) * 100; 
PctIncQdQmTot = (QTd-Om_mpg.QT) / Om_mpg.QT * 100; PctIncPVQdQmTot = (PVQTd-Om_mpg.PVQT) / Om_mpg.PVQT * 100
PctIncQdQiTot = (QTd-Oi_mpg.QT) / Oi_mpg.QT * 100; PctIncPVQdQiTot = (PVQTd-Oi_mpg.PVQT) / Oi_mpg.PVQT * 100
PctIncQmQiTot = (Om_mpg.QT-Oi_mpg.QT) / Oi_mpg.QT * 100; PctIncPVQmQiTot = (Om_mpg.PVQT-Oi_mpg.PVQT) / Oi_mpg.PVQT * 100
Op_mpg = (PctIncDdDiTot = PctIncDdDiTot, PctIncDdDmTot = PctIncDdDmTot, 
        PctIncQdQmTot = PctIncQdQmTot, PctIncPVQdQmTot = PctIncPVQdQmTot, 
        PctIncQdQiTot = PctIncQdQiTot, PctIncPVQdQiTot = PctIncPVQdQiTot,
        PctIncQmQiTot = PctIncQmQiTot, PctIncPVQmQiTot = PctIncPVQmQiTot)


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altmp.jld2" 

