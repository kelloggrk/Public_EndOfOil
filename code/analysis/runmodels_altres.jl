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
results_altres = Dict()
Threads.@threads for spec in ["p_lowres", "p_highres"]
    p = data[spec]
    results_altres[spec] = runAKSmodel(p)
end
Oi_lowres, Od_lowres, Om_lowres, Op_lowres = results_altres["p_lowres"]
Oi_highres, Od_highres, Om_highres, Op_highres = results_altres["p_highres"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altres.jld2" 

