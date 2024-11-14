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
results_altdall1 = Dict()
Threads.@threads for spec in ["p_dall03", "p_dall04", "p_dall05", "p_dall06", "p_dall07", "p_dall08", "p_dall09"]
    p = data[spec]
    results_altdall1[spec] = runAKSmodel(p)
end
Oi_dall03, Od_dall03, Om_dall03, Op_dall03 = results_altdall1["p_dall03"]
Oi_dall04, Od_dall04, Om_dall04, Op_dall04 = results_altdall1["p_dall04"]
Oi_dall05, Od_dall05, Om_dall05, Op_dall05 = results_altdall1["p_dall05"]
Oi_dall06, Od_dall06, Om_dall06, Op_dall06 = results_altdall1["p_dall06"]
Oi_dall07, Od_dall07, Om_dall07, Op_dall07 = results_altdall1["p_dall07"]
Oi_dall08, Od_dall08, Om_dall08, Op_dall08 = results_altdall1["p_dall08"]
Oi_dall09, Od_dall09, Om_dall09, Op_dall09 = results_altdall1["p_dall09"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altdall1.jld2" 

