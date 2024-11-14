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
results_n = Dict()
Threads.@threads for spec in ["p_n", "p_n_lowres", "p_n_dall03"]
    p = data[spec]
    results_n[spec] = runnoonAKSmodel(p)
end
Oi_n, Od_n, Om_n, Op_n = results_n["p_n"]
Oi_n_lowres, Od_n_lowres, Om_n_lowres, Op_n_lowres = results_n["p_n_lowres"]
Oi_n_dall03, Od_n_dall03, Om_n_dall03, Op_n_dall03 = results_n["p_n_dall03"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_n.jld2" 

