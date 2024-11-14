# ==============================================================
# Script to output all results from the models
# Outputs figures, tables, and single number tex files
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
include(dircurrent * "PlotFunctions.jl")



# ==============================================================
# Load results and input parameters
# ==============================================================
resfile_altmp = "EOOmodels_altmp.jld2"
resfile_altdec = "EOOmodels_altdec.jld2"
resfile_n = "EOOmodels_n.jld2"
resfile_altres = "EOOmodels_altres.jld2"
resfile_altdall1 = "EOOmodels_altdall1.jld2"
resfile_altdall2 = "EOOmodels_altdall2.jld2"
resfile_altdemdec = "EOOmodels_altdemdec.jld2"
resfile_altelast = "EOOmodels_altelast.jld2"

# Reference case
variable_names = ["Oi_ref", "Od_ref", "Om_ref", "Op_ref"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altmp, \"$(var_name)\")"))
end
# Non-AKS model
variable_names = ["Oi_n", "Od_n", "Om_n", "Op_n"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_n, \"$(var_name)\")"))
end
# Alternative decline rates
variable_names = ["Oi_dec30", "Od_dec30", "Om_dec30", "Op_dec30"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdec, \"$(var_name)\")"))
end
variable_names = ["Oi_dec08", "Od_dec08", "Om_dec08", "Op_dec08"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdec, \"$(var_name)\")"))
end
variable_names = ["Oi_dec06", "Od_dec06", "Om_dec06", "Op_dec06"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdec, \"$(var_name)\")"))
end
# Alternative reserves
variable_names = ["Oi_lowres", "Od_lowres", "Om_lowres", "Op_lowres"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altres, \"$(var_name)\")"))
end
variable_names = ["Oi_highres", "Od_highres", "Om_highres", "Op_highres"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altres, \"$(var_name)\")"))
end
# Alternative discount rates
variable_names = ["Oi_dall03", "Od_dall03", "Om_dall03", "Op_dall03"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall04", "Od_dall04", "Om_dall04", "Op_dall04"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall05", "Od_dall05", "Om_dall05", "Op_dall05"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall06", "Od_dall06", "Om_dall06", "Op_dall06"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall07", "Od_dall07", "Om_dall07", "Op_dall07"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall08", "Od_dall08", "Om_dall08", "Op_dall08"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall09", "Od_dall09", "Om_dall09", "Op_dall09"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall1, \"$(var_name)\")"))
end
variable_names = ["Oi_dall10", "Od_dall10", "Om_dall10", "Op_dall10"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
variable_names = ["Oi_dall12", "Od_dall12", "Om_dall12", "Op_dall12"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
variable_names = ["Oi_dall14", "Od_dall14", "Om_dall14", "Op_dall14"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
variable_names = ["Oi_dall16", "Od_dall16", "Om_dall16", "Op_dall16"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
variable_names = ["Oi_dall18", "Od_dall18", "Om_dall18", "Op_dall18"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
variable_names = ["Oi_dall20", "Od_dall20", "Om_dall20", "Op_dall20"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdall2, \"$(var_name)\")"))
end
# non-AKS model with alternative reserves and discount rates
variable_names = ["Oi_n_lowres", "Od_n_lowres", "Om_n_lowres", "Op_n_lowres"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_n, \"$(var_name)\")"))
end
variable_names = ["Oi_n_dall03", "Od_n_dall03", "Om_n_dall03", "Op_n_dall03"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_n, \"$(var_name)\")"))
end
# Faster and slower demand decline
variable_names = ["Oi_T50", "Od_T50", "Om_T50", "Op_T50"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdemdec, \"$(var_name)\")"))
end
variable_names = ["Oi_T100", "Od_T100", "Om_T100", "Op_T100"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdemdec, \"$(var_name)\")"))
end
# Delayed demand decline
variable_names = ["Oi_del05", "Od_del05", "Om_del05", "Op_del05"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdemdec, \"$(var_name)\")"))
end
# Delayed demand decline, 9% discounting
variable_names = ["Oi_del05_dall09", "Od_del05_dall09", "Om_del05_dall09", "Op_del05_dall09"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdemdec, \"$(var_name)\")"))
end
# Incomplete demand decline 
variable_names = ["Oi_dem85", "Od_dem85", "Om_dem85", "Op_dem85"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altdemdec, \"$(var_name)\")"))
end
# Russia in core OPEC
variable_names = ["Oi_R", "Od_R", "Om_R", "Op_R"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altmp, \"$(var_name)\")"))
end
# No market power
variable_names = ["Oi_nmp", "Od_nmp", "Om_nmp", "Op_nmp"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altmp, \"$(var_name)\")"))
end
# Market power goes away during decline
variable_names = ["Oi_mpg", "Od_mpg", "Om_mpg", "Op_mpg"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altmp, \"$(var_name)\")"))
end
# High and low demand elasticities
variable_names = ["Oi_delhigh", "Od_delhigh", "Om_delhigh", "Op_delhigh"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end
variable_names = ["Oi_dellow", "Od_dellow", "Om_dellow", "Op_dellow"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end
# No demand growth
variable_names = ["Oi_dng", "Od_dng", "Om_dng", "Op_dng"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end
# High demand growth
variable_names = ["Oi_dhg", "Od_dhg", "Om_dhg", "Op_dhg"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end
# Higher cost function intercepts
variable_names = ["Oi_higha", "Od_higha", "Om_higha", "Op_higha"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end
variable_names = ["Oi_highera", "Od_highera", "Om_highera", "Op_highera"]
for var_name in variable_names
    eval(Meta.parse("$(var_name) = load(dirint*resfile_altelast, \"$(var_name)\")"))
end

# Load reference case parameters
p_ref = load(dirint * "cal.jld2", "p_ref")
prefT = p_ref.T
# Load parameters for some alternative cases (for single-number tex files)
p_lowres = load(dirint * "cal.jld2", "p_lowres")
p_highres = load(dirint * "cal.jld2", "p_highres")
p_T50 = load(dirint * "cal.jld2", "p_T50")
p_T100 = load(dirint * "cal.jld2", "p_T100")
p_del05 = load(dirint * "cal.jld2", "p_del05")
p_dem85 = load(dirint * "cal.jld2", "p_dem85")
p_dec06 = load(dirint * "cal.jld2", "p_dec06")
p_R = load(dirint * "cal.jld2", "p_R")
p_delhigh = load(dirint * "cal.jld2", "p_delhigh")
p_dellow = load(dirint * "cal.jld2", "p_dellow")
p_dng = load(dirint * "cal.jld2", "p_dng")
p_dhg = load(dirint * "cal.jld2", "p_dhg")
p_highera = load(dirint * "cal.jld2", "p_highera")



# ==============================================================
# Plot aggregate results for ref case
# ==============================================================
# Paper
comboplot_ref = comboplot(Oi_ref, Od_ref, Om_ref, prefT, true)
savefig(comboplot_ref, dirfig*"comboplot_ref.pdf")
# Slides 
comboplot_ref_slides = comboplot(Oi_ref, Od_ref, Om_ref, prefT, false)
savefig(comboplot_ref_slides, dirfigs*"comboplot_ref.pdf")
comboplot_ref_slides2 = comboplot(Oi_ref, Od_ref, Om_ref, prefT, false, 2)
savefig(comboplot_ref_slides2, dirfigs*"comboplot_ref2.pdf")
comboplot_ref_slides1 = comboplot(Oi_ref, Od_ref, Om_ref, prefT, false, 1)
savefig(comboplot_ref_slides1, dirfigs*"comboplot_ref1.pdf")



# ==============================================================
# Plot region results for ref case
# ==============================================================
N = length(p_ref.q0)
legvec = [true, false, false, false]
yaxtvec = [true, false, true, false]
for region in 1:N
    rplot_ref = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT, 
        region, true, legvec[region], yaxtvec[region])
    savefig(rplot_ref, dirfig * "rplot_ref_region$region.pdf")
end
# Version for slides 
rplot_ref_slides1 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT, 
    1, false, legvec[1], yaxtvec[1])
rplot_ref_slides2 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    2, false, legvec[2], yaxtvec[2])
rplot_ref_slides3 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    3, false, legvec[3], yaxtvec[3])
rplot_ref_slides4 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    4, false, legvec[4], yaxtvec[4])
title!(rplot_ref_slides1, "Core OPEC")
title!(rplot_ref_slides2, "Non-core OPEC+")
title!(rplot_ref_slides3, "Conventional non-OPEC")
title!(rplot_ref_slides4, "Shale oil")
combined_region_ref_slides = plot(rplot_ref_slides1, rplot_ref_slides2, rplot_ref_slides3, rplot_ref_slides4, 
    layout=@layout([a b; c d]), size=(1000, 900))
savefig(combined_region_ref_slides, dirfigs*"combined_region_ref_slides.pdf")
# Version for slides, baseline only
rplot_ref_slides1 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT, 
    1, false, legvec[1], yaxtvec[1], 1)
rplot_ref_slides2 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    2, false, legvec[2], yaxtvec[2], 1)
rplot_ref_slides3 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    3, false, legvec[3], yaxtvec[3], 1)
rplot_ref_slides4 = regionplot(Oi_ref.qvec, Od_ref.qvec, Om_ref.qvec, prefT,
    4, false, legvec[4], yaxtvec[4], 1)
title!(rplot_ref_slides1, "Core OPEC")
title!(rplot_ref_slides2, "Non-core OPEC+")
title!(rplot_ref_slides3, "Conventional non-OPEC")
title!(rplot_ref_slides4, "Shale oil")
combined_region_ref_slides1 = plot(rplot_ref_slides1, rplot_ref_slides2, rplot_ref_slides3, rplot_ref_slides4, 
    layout=@layout([a b; c d]), size=(1000, 900))
savefig(combined_region_ref_slides1, dirfigs*"combined_region_ref_slides1.pdf")

# Region plot for slides, 3% discounting with no investment
rplot_n_dall03_slides1 = regionplot(Oi_n_dall03.qvec, Od_n_dall03.qvec, Om_n_dall03.qvec, prefT, 
    1, false, legvec[1], yaxtvec[1])
rplot_n_dall03_slides2 = regionplot(Oi_n_dall03.qvec, Od_n_dall03.qvec, Om_n_dall03.qvec, prefT,
    2, false, legvec[2], yaxtvec[2])
rplot_n_dall03_slides3 = regionplot(Oi_n_dall03.qvec, Od_n_dall03.qvec, Om_n_dall03.qvec, prefT,
    3, false, legvec[3], yaxtvec[3])
rplot_n_dall03_slides4 = regionplot(Oi_n_dall03.qvec, Od_n_dall03.qvec, Om_n_dall03.qvec, prefT,
    4, false, legvec[4], yaxtvec[4])
title!(rplot_n_dall03_slides1, "Core OPEC")
title!(rplot_n_dall03_slides2, "Non-core OPEC+")
title!(rplot_n_dall03_slides3, "Conventional non-OPEC")
title!(rplot_n_dall03_slides4, "Shale oil")
combined_region_n_dall03_slides = plot(rplot_n_dall03_slides1, rplot_n_dall03_slides2, rplot_n_dall03_slides3, rplot_n_dall03_slides4, 
    layout=@layout([a b; c d]), size=(1000, 900))
savefig(combined_region_n_dall03_slides, dirfigs*"combined_region_n_dall03_slides.pdf")



# ==============================================================
# Plot aggregate results for alternative specs
# ==============================================================
# Define the specs
specs = ["dec06", "dec08", "dec30", "n", "lowres", "highres", "dall09", "dall03", "n_lowres", "n_dall03", 
    "T50", "T100", "del05", "del05_dall09", "dem85"]
# Loop over the specs
for spec in specs
    # Construct variable names
    Oi = eval(Symbol("Oi_$spec"))
    Od = eval(Symbol("Od_$spec"))
    Om = eval(Symbol("Om_$spec"))    
    # Generate and save the plot
    comboplot_result = comboplot(Oi, Od, Om, prefT, true)
    savefig(comboplot_result, dirfig * "comboplot_$spec.pdf")
end



# ==============================================================
# Plot aggregate results for alternative specs, for slides
# ==============================================================
comboplot_del05_slides = comboplot(Oi_del05, Od_del05, Om_del05, prefT, false)
savefig(comboplot_del05_slides, dirfigs*"comboplot_del05.pdf")
comboplot_del05_dall09_slides = comboplot(Oi_del05_dall09, Od_del05_dall09, Om_del05_dall09, prefT, false)
savefig(comboplot_del05_dall09_slides, dirfigs*"comboplot_del05_dall09.pdf")
# Zoomed in price plots
pricezoomplot_del05_slides = priceplotzoomed(Oi_del05, Od_del05, Om_del05, prefT, false)
savefig(pricezoomplot_del05_slides, dirfigs*"pricezoomplot_del05.pdf")
pricezoomplot_del05_dall09_slides = priceplotzoomed(Oi_del05_dall09, Od_del05_dall09, Om_del05_dall09, prefT, false)
savefig(pricezoomplot_del05_dall09_slides, dirfigs*"pricezoomplot_del05_dall09.pdf")



# ==============================================================
# Plot net anticipation effect vs discount rate
# ==============================================================
xdisc = [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]
ydisc = [Op_dall03.PctIncQdQmTot, Op_dall04.PctIncQdQmTot, Op_dall05.PctIncQdQmTot, Op_dall06.PctIncQdQmTot, 
    Op_dall07.PctIncQdQmTot, Op_dall08.PctIncQdQmTot, Op_dall09.PctIncQdQmTot, Op_dall10.PctIncQdQmTot, 
    Op_dall12.PctIncQdQmTot, Op_dall14.PctIncQdQmTot, Op_dall16.PctIncQdQmTot, Op_dall18.PctIncQdQmTot, 
    Op_dall20.PctIncQdQmTot]
fsize = 12
plot(xdisc, ydisc, linewidth=2, linecolor=:black, legend=false)
xticks!(0:2:20)
xlabel!("Discount rate, %")
ylabel!("Net anticipation effect, %")
plot!(tickfontsize=fsize, titlefontsize=fsize, xguidefontsize=fsize, 
    yguidefontsize=fsize, fontfamily="Computer Modern")
savefig(dirfig * "discplot.pdf")


# ==============================================================
# Table of aggregate quantity changes - ref case and main alternatives
# ==============================================================
TabAgg  = "\\begin{tabular} {l l c c c c c c c c c c} \\midrule \\midrule \n"
TabAgg *= " & & \\textbf{(1)} & \\textbf{(2)} & \\textbf{(3)} & \\textbf{(4)} & \\textbf{(5)} & \\textbf{(6)} & \\textbf{(7)} & \\textbf{(8)} & \\textbf{(9)} & \\textbf{(10)} \\\\ \n"
TabAgg *= " & & & & \\multicolumn{2}{c}{\\textbf{Anticipated}} & \\multicolumn{2}{c}{\\textbf{Unanticipated}} & & & & \\\\ \n"
TabAgg *= " & & \\multicolumn{2}{c}{\\textbf{Baseline}} & \\multicolumn{2}{c}{\\textbf{decline}} & \\multicolumn{2}{c}{\\textbf{decline}} &
    \\underline{\\textbf{(3)-(1)}} & \\underline{\\textbf{(4)-(2)}} & \\underline{\\textbf{(3)-(5)}} & \\underline{\\textbf{(4)-(6)}} \\\\ \n"
TabAgg *= " & \\multicolumn{1}{c}{\\textbf{Model}} & \\textbf{Q} & \\textbf{PVQ} & \\textbf{Q} & \\textbf{PVQ} & \\textbf{Q} & \\textbf{PVQ} & 
    \\textbf{(1)} & \\textbf{(2)} & \\textbf{(5)} & \\textbf{(6)} \\\\ \n"
TabAgg *= "\\midrule \n"
TabAgg *= @sprintf("\\textbf{1.} & \\textbf{Reference case} & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_ref.Q), sum(Oi_ref.PVQ), sum(Od_ref.Q), sum(Od_ref.PVQ), sum(Om_ref.Q), sum(Om_ref.PVQ),
    Op_ref.PctIncQdQiTot, Op_ref.PctIncPVQdQiTot, Op_ref.PctIncQdQmTot, Op_ref.PctIncPVQdQmTot)
TabAgg *= "\\midrule \n"
TabAgg *= @sprintf("\\multicolumn{4}{l}{\\textbf{Alternative declines, reserves, and discount rates}} & & & & & & & & \\\\ \n")
TabAgg *= @sprintf("2. & All regions decline at 8\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dec08.Q), sum(Oi_dec08.PVQ), sum(Od_dec08.Q), sum(Od_dec08.PVQ), sum(Om_dec08.Q), sum(Om_dec08.PVQ),
    Op_dec08.PctIncQdQiTot, Op_dec08.PctIncPVQdQiTot, Op_dec08.PctIncQdQmTot, Op_dec08.PctIncPVQdQmTot)
TabAgg *= @sprintf("3. & All regions decline at 6\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dec08.Q), sum(Oi_dec06.PVQ), sum(Od_dec06.Q), sum(Od_dec06.PVQ), sum(Om_dec06.Q), sum(Om_dec06.PVQ),
    Op_dec06.PctIncQdQiTot, Op_dec06.PctIncPVQdQiTot, Op_dec06.PctIncQdQmTot, Op_dec06.PctIncPVQdQmTot)
TabAgg *= @sprintf("4. & All regions decline at 30\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dec30.Q), sum(Oi_dec30.PVQ), sum(Od_dec30.Q), sum(Od_dec30.PVQ), sum(Om_dec30.Q), sum(Om_dec30.PVQ),
    Op_dec30.PctIncQdQiTot, Op_dec30.PctIncPVQdQiTot, Op_dec30.PctIncQdQmTot, Op_dec30.PctIncPVQdQmTot)
TabAgg *= @sprintf("5. & No investment & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_n.Q), sum(Oi_n.PVQ), sum(Od_n.Q), sum(Od_n.PVQ), sum(Om_n.Q), sum(Om_n.PVQ),
    Op_n.PctIncQdQiTot, Op_n.PctIncPVQdQiTot, Op_n.PctIncQdQmTot, Op_n.PctIncPVQdQmTot)
TabAgg *= @sprintf("6. & High reserves (%.0f billion bbl) & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(p_highres.x0), sum(Oi_highres.Q), sum(Oi_highres.PVQ), sum(Od_highres.Q), sum(Od_highres.PVQ), sum(Om_highres.Q), sum(Om_highres.PVQ),
    Op_highres.PctIncQdQiTot, Op_highres.PctIncPVQdQiTot, Op_highres.PctIncQdQmTot, Op_highres.PctIncPVQdQmTot)
TabAgg *= @sprintf("7. & Low reserves (%.0f billion bbl) & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(p_lowres.x0), sum(Oi_lowres.Q), sum(Oi_lowres.PVQ), sum(Od_lowres.Q), sum(Od_lowres.PVQ), sum(Om_lowres.Q), sum(Om_lowres.PVQ),
    Op_lowres.PctIncQdQiTot, Op_lowres.PctIncPVQdQiTot, Op_lowres.PctIncQdQmTot, Op_lowres.PctIncPVQdQmTot)
TabAgg *= @sprintf("8. & All resource types discount at 9\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dall09.Q), sum(Oi_dall09.PVQ), sum(Od_dall09.Q), sum(Od_dall09.PVQ), sum(Om_dall09.Q), sum(Om_dall09.PVQ),
    Op_dall09.PctIncQdQiTot, Op_dall09.PctIncPVQdQiTot, Op_dall09.PctIncQdQmTot, Op_dall09.PctIncPVQdQmTot)
TabAgg *= @sprintf("9. & All resource types discount at 3\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dall03.Q), sum(Oi_dall03.PVQ), sum(Od_dall03.Q), sum(Od_dall03.PVQ), sum(Om_dall03.Q), sum(Om_dall03.PVQ),
    Op_dall03.PctIncQdQiTot, Op_dall03.PctIncPVQdQiTot, Op_dall03.PctIncQdQmTot, Op_dall03.PctIncPVQdQmTot)
TabAgg *= @sprintf("10. & No investment + low reserves & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_n_lowres.Q), sum(Oi_n_lowres.PVQ), sum(Od_n_lowres.Q), sum(Od_n_lowres.PVQ), sum(Om_n_lowres.Q), sum(Om_n_lowres.PVQ),
    Op_n_lowres.PctIncQdQiTot, Op_n_lowres.PctIncPVQdQiTot, Op_n_lowres.PctIncQdQmTot, Op_n_lowres.PctIncPVQdQmTot)
TabAgg *= @sprintf("11. & No investment + discounting at 3\\%% & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_n_dall03.Q), sum(Oi_n_dall03.PVQ), sum(Od_n_dall03.Q), sum(Od_n_dall03.PVQ), sum(Om_n_dall03.Q), sum(Om_n_dall03.PVQ),
    Op_n_dall03.PctIncQdQiTot, Op_n_dall03.PctIncPVQdQiTot, Op_n_dall03.PctIncQdQmTot, Op_n_dall03.PctIncPVQdQmTot)
TabAgg *= "\\midrule \n"
TabAgg *= @sprintf("\\multicolumn{4}{l}{\\textbf{Alternative demand declines}} & & & & & & & & \\\\ \n")
TabAgg *= @sprintf("12. & 50-year demand decline & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_T50.Q), sum(Oi_T50.PVQ), sum(Od_T50.Q), sum(Od_T50.PVQ), sum(Om_T50.Q), sum(Om_T50.PVQ),
    Op_T50.PctIncQdQiTot, Op_T50.PctIncPVQdQiTot, Op_T50.PctIncQdQmTot, Op_T50.PctIncPVQdQmTot)
TabAgg *= @sprintf("13. & 100-year demand decline & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_T100.Q), sum(Oi_T100.PVQ), sum(Od_T100.Q), sum(Od_T100.PVQ), sum(Om_T100.Q), sum(Om_T100.PVQ),
    Op_T100.PctIncQdQiTot, Op_T100.PctIncPVQdQiTot, Op_T100.PctIncQdQmTot, Op_T100.PctIncPVQdQmTot)
TabAgg *= @sprintf("14. & Decline delayed 5 years & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_del05.Q), sum(Oi_del05.PVQ), sum(Od_del05.Q), sum(Od_del05.PVQ), sum(Om_del05.Q), sum(Om_del05.PVQ),
    Op_del05.PctIncQdQiTot, Op_del05.PctIncPVQdQiTot, Op_del05.PctIncQdQmTot, Op_del05.PctIncPVQdQmTot)
TabAgg *= @sprintf("15. & Decline delayed 5 years; 9\\%% discounting & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_del05_dall09.Q), sum(Oi_del05_dall09.PVQ), sum(Od_del05_dall09.Q), sum(Od_del05_dall09.PVQ), sum(Om_del05_dall09.Q), sum(Om_del05_dall09.PVQ),
    Op_del05_dall09.PctIncQdQiTot, Op_del05_dall09.PctIncPVQdQiTot, Op_del05_dall09.PctIncQdQmTot, Op_del05_dall09.PctIncPVQdQmTot)
TabAgg *= @sprintf("16. & 15\\%% of demand remains after decline & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dem85.Q), sum(Oi_dem85.PVQ), sum(Od_dem85.Q), sum(Od_dem85.PVQ), sum(Om_dem85.Q), sum(Om_dem85.PVQ),
    Op_dem85.PctIncQdQiTot, Op_dem85.PctIncPVQdQiTot, Op_dem85.PctIncQdQmTot, Op_dem85.PctIncPVQdQmTot)
TabAgg *= "\\midrule \n"
TabAgg *= "\\end{tabular}"
write(dirtab * "TabAgg.tex", TabAgg)



# ==============================================================
# Table of aggregate quantity changes - other alternatives
# ==============================================================
TabAgg  = "\\begin{tabular} {l l c c c c c c c c c c} \\midrule \\midrule \n"
TabAgg *= " & & \\textbf{(1)} & \\textbf{(2)} & \\textbf{(3)} & \\textbf{(4)} & \\textbf{(5)} & \\textbf{(6)} & \\textbf{(7)} & \\textbf{(8)} & \\textbf{(9)} & \\textbf{(10)} \\\\ \n"
TabAgg *= " & & & & \\multicolumn{2}{c}{\\textbf{Anticipated}} & \\multicolumn{2}{c}{\\textbf{Unanticipated}} & & & & \\\\ \n"
TabAgg *= " & & \\multicolumn{2}{c}{\\textbf{Baseline}} & \\multicolumn{2}{c}{\\textbf{decline}} & \\multicolumn{2}{c}{\\textbf{decline}} &
    \\underline{\\textbf{(3)-(1)}} & \\underline{\\textbf{(4)-(2)}} & \\underline{\\textbf{(3)-(5)}} & \\underline{\\textbf{(4)-(6)}} \\\\ \n"
TabAgg *= " & \\multicolumn{1}{c}{\\textbf{Model}} & \\textbf{Q} & \\textbf{PVQ} & \\textbf{Q} & \\textbf{PVQ} & \\textbf{Q} & \\textbf{PVQ} & 
    \\textbf{(1)} & \\textbf{(2)} & \\textbf{(5)} & \\textbf{(6)} \\\\ \n"
TabAgg *= "\\midrule \n"
TabAgg *= @sprintf("\\textbf{1.} & \\textbf{Reference case} & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_ref.Q), sum(Oi_ref.PVQ), sum(Od_ref.Q), sum(Od_ref.PVQ), sum(Om_ref.Q), sum(Om_ref.PVQ),
    Op_ref.PctIncQdQiTot, Op_ref.PctIncPVQdQiTot, Op_ref.PctIncQdQmTot, Op_ref.PctIncPVQdQmTot)
TabAgg *= "\\midrule \n"
TabAgg *= @sprintf("\\multicolumn{4}{l}{\\textbf{Alternative specifications}} & & & & & & & & \\\\ \n")
TabAgg *= @sprintf("17. & Russia in core OPEC, 9\\%% discounting & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_R.Q), sum(Oi_R.PVQ), sum(Od_R.Q), sum(Od_R.PVQ), sum(Om_R.Q), sum(Om_R.PVQ),
    Op_R.PctIncQdQiTot, Op_R.PctIncPVQdQiTot, Op_R.PctIncQdQmTot, Op_R.PctIncPVQdQmTot)
TabAgg *= @sprintf("18. & No market power & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_nmp.Q), sum(Oi_nmp.PVQ), sum(Od_nmp.Q), sum(Od_nmp.PVQ), sum(Om_nmp.Q), sum(Om_nmp.PVQ),
    Op_nmp.PctIncQdQiTot, Op_nmp.PctIncPVQdQiTot, Op_nmp.PctIncQdQmTot, Op_nmp.PctIncPVQdQmTot)
TabAgg *= @sprintf("19. & No market power in anticipated decline & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_mpg.Q), sum(Oi_mpg.PVQ), sum(Od_mpg.Q), sum(Od_mpg.PVQ), sum(Om_mpg.Q), sum(Om_mpg.PVQ),
    Op_mpg.PctIncQdQiTot, Op_mpg.PctIncPVQdQiTot, Op_mpg.PctIncQdQmTot, Op_mpg.PctIncPVQdQmTot)
TabAgg *= @sprintf("20. & High demand elasticity (%.1f) & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    p_delhigh.delast, sum(Oi_delhigh.Q), sum(Oi_delhigh.PVQ), sum(Od_delhigh.Q), sum(Od_delhigh.PVQ), sum(Om_delhigh.Q), sum(Om_delhigh.PVQ),
    Op_delhigh.PctIncQdQiTot, Op_delhigh.PctIncPVQdQiTot, Op_delhigh.PctIncQdQmTot, Op_delhigh.PctIncPVQdQmTot)
TabAgg *= @sprintf("21. & Low demand elasticity (%.1f) & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    p_dellow.delast, sum(Oi_dellow.Q), sum(Oi_dellow.PVQ), sum(Od_dellow.Q), sum(Od_dellow.PVQ), sum(Om_dellow.Q), sum(Om_dellow.PVQ),
    Op_dellow.PctIncQdQiTot, Op_dellow.PctIncPVQdQiTot, Op_dellow.PctIncQdQmTot, Op_dellow.PctIncPVQdQmTot)
TabAgg *= @sprintf("22. & No demand growth & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_dng.Q), sum(Oi_dng.PVQ), sum(Od_dng.Q), sum(Od_dng.PVQ), sum(Om_dng.Q), sum(Om_dng.PVQ),
    Op_dng.PctIncQdQiTot, Op_dng.PctIncPVQdQiTot, Op_dng.PctIncQdQmTot, Op_dng.PctIncPVQdQmTot)
TabAgg *= @sprintf("23. & Higher demand growth (%.1f\\%% per year) & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    p_dhg.dgr_ann*100, sum(Oi_dhg.Q), sum(Oi_dhg.PVQ), sum(Od_dhg.Q), sum(Od_dhg.PVQ), sum(Om_dhg.Q), sum(Om_dhg.PVQ),
    Op_dhg.PctIncQdQiTot, Op_dhg.PctIncPVQdQiTot, Op_dhg.PctIncQdQmTot, Op_dhg.PctIncPVQdQmTot)
TabAgg *= @sprintf("24. & High cost function intercepts & %d & %d & %d & %d & %d & %d & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n", 
    sum(Oi_highera.Q), sum(Oi_highera.PVQ), sum(Od_highera.Q), sum(Od_highera.PVQ), sum(Om_highera.Q), sum(Om_highera.PVQ),
    Op_highera.PctIncQdQiTot, Op_highera.PctIncPVQdQiTot, Op_highera.PctIncQdQmTot, Op_highera.PctIncPVQdQmTot)
TabAgg *= "\\midrule \n"
TabAgg *= "\\end{tabular}"
write(dirtab * "TabAgg_alt.tex", TabAgg)



# ==============================================================
# Output single-number tex files for reference case
# ==============================================================
# Last drilling year for non core OPEC at baseline
Ys = 2023
open(dirsnt*"sim/Te_noncore.tex", "w") do file
    write(file, string(Int(floor(Ys+maximum(Oi_ref.Te[2:N])/prefT))))
end
# Last drilling year for core OPEC at baseline
open(dirsnt*"sim/Te_core.tex", "w") do file
    write(file, string(Int(floor(Ys+maximum(Oi_ref.Te[1:N])/prefT))))
end
# Scarcity rents in $/bbl
for region in 1:4
    open(dirsnt * "sim/mu0_ref_$(region).tex", "w") do file
        write(file, @sprintf("%.2f",Oi_ref.mu0_bbl[region]))
    end
end
# Max non-core OPEC scarcity rent
open(dirsnt*"sim/mu0_ref_noncore.tex", "w") do file
    write(file, string(Int(ceil(maximum(Oi_ref.mu0_bbl[2:N]), digits=0))))
end
# Years of demand decline
open(dirsnt*"sim/Ttz_ref.tex", "w") do file
    write(file, string(Int(round(p_ref.timetozero, digits=0))))
end
# Total extraction under anticpated decline
open(dirsnt * "sim/Qd_ref_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_ref.Q), digits=0))))
end
# Initial investment at baseline and under the anticipated declinhe
open(dirsnt * "sim/di1_ref.tex", "w") do file
    write(file, @sprintf("%.2f", sum(Oi_ref.dvec[1,:])))
end
open(dirsnt * "sim/dd1_ref.tex", "w") do file
    write(file, @sprintf("%.2f", sum(Od_ref.dvec[1,:])))
end
# Percent changes in total
open(dirsnt * "sim/PctDecDdDi_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncDdDiTot))
end
open(dirsnt * "sim/PctDecQdQi_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncQdQiTot))
end
open(dirsnt * "sim/PctDecDdDm_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncDdDmTot))
end
open(dirsnt * "sim/PctDecQdQm_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecPVQdQi_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncPVQdQiTot))
end
open(dirsnt * "sim/PctDecPVQdQm_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_ref.PctIncPVQdQmTot))
end
# Percent changes in first period drilling
PctIncdd0dm0_ref_tot = (sum(vec(Od_ref.dvec[1,:]))-sum(vec(Om_ref.dvec[1,:])))/sum(vec(Om_ref.dvec[1,:]))*100
open(dirsnt * "sim/PctDecdd0dm0_ref_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -PctIncdd0dm0_ref_tot))
end
# Percent changes in production by region
for region in 1:4
    open(dirsnt * "sim/PctDecQdQi_ref_$(region).tex", "w") do file
        write(file, @sprintf("%.1f", -(Od_ref.Q[region]-Oi_ref.Q[region])
            /Oi_ref.Q[region]*100))
    end
end



# ==============================================================
# Output single-number tex files for main alternative specs
# ==============================================================
# Percent changes in total
open(dirsnt * "sim/PctDecQdQm_dec06_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dec06.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dec08_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dec08.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dec30_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dec30.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_n_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_n.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_highres_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_highres.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_lowres_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_lowres.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dall09_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dall09.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dall03_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_dall03.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_n_lowres_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_n_lowres.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_n_dall03_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_n_dall03.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_T50_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_T50.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_T100_tot.tex", "w") do file
    write(file, @sprintf("%.1f", Op_T100.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_del05_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_del05.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_del05_dall09_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_del05_dall09.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dem85_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dem85.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecPVQdQm_dem85_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dem85.PctIncPVQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_R_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_R.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_nmp_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_nmp.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_mpg_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_mpg.PctIncQdQmTot))
end
open(dirsnt * "sim/PctDecQdQm_dng_tot.tex", "w") do file
    write(file, @sprintf("%.1f", -Op_dng.PctIncQdQmTot))
end
# Scarcity rents
open(dirsnt * "sim/mu0_lowres_3.tex", "w") do file
    write(file, @sprintf("%.2f", Oi_lowres.mu0_bbl[3]))
end
open(dirsnt * "sim/mu0_dall03_3.tex", "w") do file
    write(file, @sprintf("%.2f", Oi_dall03.mu0_bbl[3]))
end
for region in 1:4
    open(dirsnt * "sim/mu0_n_lowres_$(region).tex", "w") do file
        write(file, @sprintf("%.2f",Oi_n_lowres.mu0[region]))
    end
end
for region in 1:4
    open(dirsnt * "sim/mu0_n_dall03_$(region).tex", "w") do file
        write(file, @sprintf("%.2f",Oi_n_dall03.mu0[region]))
    end
end
# Discount rate at which net anticipation is lowest 
ymin = minimum(ydisc)
discind = findall(x -> x == ymin, ydisc) .+ 2
open(dirsnt * "sim/discminanticip.tex", "w") do file
    write(file, @sprintf("%.0f", discind[1]))
end
# Years of demand decline
open(dirsnt*"sim/Ttz_T50.tex", "w") do file
    write(file, string(Int(round(p_T50.timetozero, digits=0))))
end
open(dirsnt*"sim/Ttz_T100.tex", "w") do file
    write(file, string(Int(round(p_T100.timetozero, digits=0))))
end
# Delayed demand
open(dirsnt * "sim/Tdel_del05.tex", "w") do file
    write(file, string(Int(round(p_del05.shiftdelin, digits=0))))
end
open(dirsnt * "sim/Ttz_del05.tex", "w") do file
    write(file, string(Int(round(p_del05.timetozero-p_del05.shiftdelin, digits=0))))
end
# Incomplete demand decline
open(dirsnt * "sim/Drem_dem85.tex", "w") do file
    write(file, string(Int(round(p_dem85.drem*100))))
end
open(dirsnt * "sim/maxY_dem85.tex", "w") do file
    write(file, string(Int(round(p_dem85.maxY))))
end
# Total extraction
open(dirsnt * "sim/Qd_dall09_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_dall09.Q), digits=0))))
end
open(dirsnt * "sim/Qd_T50_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_T50.Q), digits=0))))
end
open(dirsnt * "sim/Qd_T100_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_T100.Q), digits=0))))
end
open(dirsnt * "sim/Qd_dem85_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_dem85.Q), digits=0))))
end
open(dirsnt * "sim/Qd_R_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_R.Q), digits=0))))
end
open(dirsnt * "sim/Qd_mpg_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_mpg.Q), digits=0))))
end
# PV extraction
open(dirsnt * "sim/PVQd_dem85_tot.tex", "w") do file
    write(file, string(Int(round(sum(Od_dem85.PVQ), digits=0))))
end
# Russian production and reserves
open(dirsnt * "sim/q0Russia.tex", "w") do file
    write(file, string(round(p_R.q0[1]-p_ref.q0[1], digits=1)))
end
open(dirsnt * "sim/x0Russia.tex", "w") do file
    write(file, string(Int(round(p_R.x0[1]-p_ref.x0[1], digits=0))))
end
# Alternative demand elasticities
open(dirsnt * "sim/delast_high.tex", "w") do file
    write(file, string(round(p_delhigh.delast, digits=1)))
end
open(dirsnt * "sim/delast_low.tex", "w") do file
    write(file, string(round(p_dellow.delast, digits=1)))
end
# Higher demand growth rate
open(dirsnt * "sim/highdemgrowth.tex", "w") do file
    write(file, string(round(p_dhg.dgr_ann*100, digits=1)))
end
# Alternative cost function intercepts
for i in 1:4
    open(dirsnt*"sim/alpha0_highera_$(i).tex", "w") do file
        write(file, string(Int(round(p_highera.alpha[i], digits=0))))
    end
end
# Oil prices in year 5 of the del05_dall09 case
open(dirsnt * "sim/p5i_del05_dall09.tex", "w") do file
    write(file, @sprintf("%.2f", Oi_del05_dall09.pvec[10]))
end
open(dirsnt * "sim/p5d_del05_dall09.tex", "w") do file
    write(file, @sprintf("%.2f", Od_del05_dall09.pvec[10]))
end

