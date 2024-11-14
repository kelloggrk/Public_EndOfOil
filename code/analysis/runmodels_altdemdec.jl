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
results_altdemdec = Dict()
Threads.@threads for spec in ["p_T50", "p_T100", "p_del05", "p_del05_dall09", "p_dem85"]
    p = data[spec]
    results_altdemdec[spec] = runAKSmodel(p)
end
Oi_T50, Od_T50, Om_T50, Op_T50 = results_altdemdec["p_T50"]
Oi_T100, Od_T100, Om_T100, Op_T100 = results_altdemdec["p_T100"]
Oi_del05, Od_del05, Om_del05, Op_del05 = results_altdemdec["p_del05"]
Oi_del05_dall09, Od_del05_dall09, Om_del05_dall09, Op_del05_dall09 = results_altdemdec["p_del05_dall09"]
Oi_dem85, Od_dem85, Om_dem85, Op_dem85 = results_altdemdec["p_dem85"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altdemdec.jld2" 

