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
results_altdall2 = Dict()
Threads.@threads for spec in ["p_dall10", "p_dall12", "p_dall14", "p_dall16", "p_dall18", "p_dall20"]
    p = data[spec]
    results_altdall2[spec] = runAKSmodel(p)
end
Oi_dall10, Od_dall10, Om_dall10, Op_dall10 = results_altdall2["p_dall10"]
Oi_dall12, Od_dall12, Om_dall12, Op_dall12 = results_altdall2["p_dall12"]
Oi_dall14, Od_dall14, Om_dall14, Op_dall14 = results_altdall2["p_dall14"]
Oi_dall16, Od_dall16, Om_dall16, Op_dall16 = results_altdall2["p_dall16"]
Oi_dall18, Od_dall18, Om_dall18, Op_dall18 = results_altdall2["p_dall18"]
Oi_dall20, Od_dall20, Om_dall20, Op_dall20 = results_altdall2["p_dall20"]



# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altdall2.jld2" 

