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
results_altelast = Dict()
Threads.@threads for spec in ["p_delhigh", "p_dellow", "p_dng", "p_dhg", "p_higha", "p_highera"]
    p = data[spec]
    results_altelast[spec] = runAKSmodel(p)
end
Oi_delhigh, Od_delhigh, Om_delhigh, Op_delhigh = results_altelast["p_delhigh"]
Oi_dellow, Od_dellow, Om_dellow, Op_dellow = results_altelast["p_dellow"]
Oi_dng, Od_dng, Om_dng, Op_dng = results_altelast["p_dng"]
Oi_dhg, Od_dhg, Om_dhg, Op_dhg = results_altelast["p_dhg"]
Oi_higha, Od_higha, Om_higha, Op_higha = results_altelast["p_higha"]
Oi_highera, Od_highera, Om_highera, Op_highera = results_altelast["p_highera"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altelast.jld2" 

