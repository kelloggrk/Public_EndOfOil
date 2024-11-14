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
data_altdec = load(dirint * "cal.jld2")


# ==============================================================
# Run reference case and alternative decline rates
# ==============================================================
results_altdec = Dict()
Threads.@threads for spec in ["p_ref", "p_dec30", "p_dec08", "p_dec06"]
    p = data_altdec[spec]
    results_altdec[spec] = runAKSmodel(p)
end
Oi_ref, Od_ref, Om_ref, Op_ref = results_altdec["p_ref"]
Oi_dec30, Od_dec30, Om_dec30, Op_dec30 = results_altdec["p_dec30"]
Oi_dec08, Od_dec08, Om_dec08, Op_dec08 = results_altdec["p_dec08"]
Oi_dec06, Od_dec06, Om_dec06, Op_dec06 = results_altdec["p_dec06"]


# ==============================================================
# Save all output
# ==============================================================
@save dirint * "EOOmodels_altdec.jld2" 

