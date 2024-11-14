# ==============================================================
# Outputs all calibrated parameters to paper as tables and
# single-number tex files
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
# Reference case
datacal = load(dirint * "cal.jld2")
p_ref0 = datacal["p_ref0"]; p_ref = datacal["p_ref"];
# Alternative cases
p_n = datacal["p_n"]
p_dec30 = datacal["p_dec30"] ; p_dec08 = datacal["p_dec08"]; p_dec06 = datacal["p_dec06"]
p_lowres = datacal["p_lowres"] ; p_highres = datacal["p_highres"]
p_dall03 = datacal["p_dall03"] ; p_dall09 = datacal["p_dall09"]
p_n_lowres = datacal["p_n_lowres"] ; p_n_dall03 = datacal["p_n_dall03"]
p_R = datacal["p_R"] ; p_nmp = datacal["p_nmp"]
p_delhigh = datacal["p_delhigh"] ; p_dellow = datacal["p_dellow"]
p_dng = datacal["p_dng"] ; p_dhg = datacal["p_dhg"]
p_higha = datacal["p_higha"] ; p_highera = datacal["p_highera"]

# elasticities
dataelast = load(dirint * "elast.jld2")
elast = dataelast["elast"]; elasttot = dataelast["elasttot"]
ce_ref = dataelast["ce_ref"]; ce_n = dataelast["ce_n"]
ce_dec30 = dataelast["ce_dec30"]; ce_dec08 = dataelast["ce_dec08"]; ce_dec06 = dataelast["ce_dec06"]
ce_lowres = dataelast["ce_lowres"]; ce_highres = dataelast["ce_highres"]
ce_dall03 = dataelast["ce_dall03"]; ce_dall09 = dataelast["ce_dall09"]
ce_n_lowres = dataelast["ce_n_lowres"]; ce_n_dall03 = dataelast["ce_n_dall03"]
ce_R = dataelast["ce_R"]; ce_nmp = dataelast["ce_nmp"]
ce_delhigh = dataelast["ce_delhigh"]; ce_dellow = dataelast["ce_dellow"]
ce_dng = dataelast["ce_dng"]; ce_dhg = dataelast["ce_dhg"]
ce_higha = dataelast["ce_higha"]; ce_highera = dataelast["ce_highera"]




# ==============================================================
# Run baseline reference case model
# ==============================================================
# Run non-AKS model to get guess of price vector and shadow value
m_temp = EOO_N(αin=p_ref.alpha, γin=p_ref.gamma, x0in=p_ref.x0, 
    delast=p_ref.delast, Pref=p_ref.Pref, dref=p_ref.dref, 
    dgr_ann=p_ref.dgr_ann, timetozero=p_ref.timetozero, 
    r_ann=p_ref.r_ann, Er_ann=p_ref.Er_ann, mp=true, T=p_ref.T)
_, pveci_temp, _, mu0i_temp, _, _ = optpath(m_temp, false, 0., x0(m_temp))
# Now run AKS model
m = EOO_AKS(eoo_n=m_temp, lambday=p_ref.lambday, q0=p_ref.q0, pguess=pveci_temp)
Di, Qi, pveci, dveci, qveci, mu0i, muveci, Tei = 
    optpath(m, q0(m), false, 1, x0(m), mu0i_temp .* RRC(m))
mc0 = mc(m, vec(dveci[1,:]))             # initial MC
mc0_bbl = mc(m, vec(dveci[1,:])) ./ RRC(m)            # initial MC in $/bbl

# Compute initial investment rates
d0 = p_ref.q0.*(dgr(m).+λ(m))



# ==============================================================
# Table of calibrated parameters
# ==============================================================
TabCal  = "\\begin{tabular} {l c c c c c c} \\midrule \\midrule \n"
TabCal *= " & &  & \\multicolumn{4}{c}{\\textbf{Value by type}} \\\\ \n"
TabCal *= "\\textbf{Parameter} &  & \\textbf{Units} & \\textbf{I} & \\textbf{II} & \\textbf{III} & \\textbf{IV} \\\\ \n"
TabCal *= "\\midrule \n"
TabCal *= "\\multicolumn{7}{l}{\\textbf{Demand parameters}} \\\\ \n"
TabCal *= @sprintf("\\hspace{4pt} Initial 2023 oil price & \$P_0\$ & \\\$/bbl & \\multicolumn{4}{c}{%.2f} \\\\ \n", p_ref.Pref)
TabCal *= @sprintf("\\hspace{4pt} Initial 2023 consumption & \$Q_0\$ & mmbbl/d & \\multicolumn{4}{c}{%.1f} \\\\ \n", sum(p_ref.q0))
TabCal *= @sprintf("\\hspace{4pt} Elasticity at \$(P_0,Q_0)\$ &  &  & \\multicolumn{4}{c}{%.1f} \\\\ \n", p_ref.delast)
TabCal *= @sprintf("\\hspace{4pt} 2023 choke price\$^*\$ &  & \\\$/bbl & \\multicolumn{4}{c}{%.2f} \\\\ \n", dintP(m))
TabCal *= @sprintf("\\hspace{4pt} Growth rate through 2030 &  & per year & \\multicolumn{4}{c}{%.4f} \\\\ \n", p_ref.dgr_ann)
TabCal *= "\\midrule \n"
TabCal *= "\\multicolumn{7}{l}{\\textbf{Supply parameters}} \\\\ \n"
TabCal *= @sprintf("\\hspace{4pt} Initial 2023 production & \$q_{i0}\$ & mmbbl/d & %.1f & %.1f & %.1f & %.1f \\\\ \n", 
    p_ref.q0[1], p_ref.q0[2], p_ref.q0[3], p_ref.q0[4])
TabCal *= @sprintf("\\hspace{4pt} Production decline rate & \$\\lambda_i\$ & per year & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_ref.lambday[1], p_ref.lambday[2], p_ref.lambday[3], p_ref.lambday[4])
TabCal *= @sprintf("\\hspace{4pt} Discount rate & \$r_i\$ & per year & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_ref.r_ann[1], p_ref.r_ann[2], p_ref.r_ann[3], p_ref.r_ann[4])
TabCal *= @sprintf("\\hspace{4pt} Reserves & \$x_{i0}\$ & billion bbl & %.0f & %.0f & %.0f & %.0f \\\\ \n", 
    p_ref.x0[1], p_ref.x0[2], p_ref.x0[3], p_ref.x0[4])
TabCal *= @sprintf("\\hspace{4pt} Initial 2023 drilling rate\$^*\$ & \$d_{i0}\$ & mmbbl/d & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    d0[1], d0[2], d0[3], d0[4])
TabCal *= @sprintf("\\hspace{4pt} Drilling cost intercept & \$a_i\$ & \\\$/bbl & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_ref.alpha[1], p_ref.alpha[2], p_ref.alpha[3], p_ref.alpha[4])
TabCal *= @sprintf("\\hspace{4pt} Drilling cost intercept\$^*\$ & \$\\alpha_i\$ & \\\$bn per mmbbl/d & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    α(m)[1], α(m)[2], α(m)[3], α(m)[4])
TabCal *= @sprintf("\\hspace{4pt} Drilling cost slope & \$g_i\$ & \\\$/bbl per mmbbl/d & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_ref.gamma[1]/λ(m)[1], p_ref.gamma[2]/λ(m)[2], p_ref.gamma[3]/λ(m)[3], p_ref.gamma[4]/λ(m)[4])
TabCal *= @sprintf("\\hspace{4pt} Drilling cost slope\$^*\$ & \$\\gamma_i\$ & \\\$bn per (mmbbl/d)\$^2\$ & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    γ(m)[1], γ(m)[2], γ(m)[3], γ(m)[4])
TabCal *= @sprintf("\\hspace{4pt} Drilling cost elasticity\$^*\$ &  &  & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_ref[1], ce_ref[2], ce_ref[3], ce_ref[4])
TabCal *= "\\midrule \n"
TabCal *= "\\end{tabular}"
write(dirtab * "TabCal.tex", TabCal)



# ==============================================================
# Table of drilling cost elasticities
# ==============================================================
TabCE  = "\\begin{tabular} {l l c c c c} \\midrule \\midrule \n"
TabCE *= " & &  \\multicolumn{4}{c}{\\textbf{Value by type}} \\\\ \n"
TabCE *= "\\multicolumn{2}{l}{\\textbf{Parameter}}  & \\textbf{I} & \\textbf{II} & \\textbf{III} & \\textbf{IV} \\\\ \n"
TabCE *= "\\midrule \n"
TabCE *= @sprintf("1. & Reference case & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_ref[1], ce_ref[2], ce_ref[3], ce_ref[4])
TabCE *= @sprintf("2. & All regions decline at 8\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dec08[1], ce_dec08[2], ce_dec08[3], ce_dec08[4])
TabCE *= @sprintf("3. & All regions decline at 6\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dec06[1], ce_dec06[2], ce_dec06[3], ce_dec06[4])
TabCE *= @sprintf("4. & All regions decline at 30\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dec30[1], ce_dec30[2], ce_dec30[3], ce_dec30[4])
TabCE *= @sprintf("5. & No investment & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_n[1], ce_n[2], ce_n[3], ce_n[4])
TabCE *= @sprintf("6. & High reserves (%.0f billion bbl) & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    sum(p_highres.x0), ce_highres[1], ce_highres[2], ce_highres[3], ce_highres[4])
TabCE *= @sprintf("7. & Low reserves (%.0f billion bbl) & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    sum(p_lowres.x0), ce_lowres[1], ce_lowres[2], ce_lowres[3], ce_lowres[4])
TabCE *= @sprintf("8. & All resource types discount at 9\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dall09[1], ce_dall09[2], ce_dall09[3], ce_dall09[4])
TabCE *= @sprintf("9. & All resource types discount at 3\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dall03[1], ce_dall03[2], ce_dall03[3], ce_dall03[4])
TabCE *= @sprintf("10. & No investment + low reserves & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_n_lowres[1], ce_n_lowres[2], ce_n_lowres[3], ce_n_lowres[4])
TabCE *= @sprintf("11. & No investment + discounting at 3\\%% & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_n_dall03[1], ce_n_dall03[2], ce_n_dall03[3], ce_n_dall03[4])
TabCE *= @sprintf("17. & Russia in core OPEC, 9\\%% discounting & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_R[1], ce_R[2], ce_R[3], ce_R[4])
TabCE *= @sprintf("18. & No market power & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_nmp[1], ce_nmp[2], ce_nmp[3], ce_nmp[4])
TabCE *= @sprintf("20. & High demand elasticity (%.1f) & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_delhigh.delast, ce_delhigh[1], ce_delhigh[2], ce_delhigh[3], ce_delhigh[4])
TabCE *= @sprintf("21. & Low demand elasticity (%.1f) & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_dellow.delast, ce_dellow[1], ce_dellow[2], ce_dellow[3], ce_dellow[4])
TabCE *= @sprintf("22. & No demand growth & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_dng[1], ce_dng[2], ce_dng[3], ce_dng[4])
TabCE *= @sprintf("23. & Higher demand growth (%.1f\\%% per year) & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    p_dhg.dgr_ann*100, ce_dhg[1], ce_dhg[2], ce_dhg[3], ce_dhg[4])
TabCE *= @sprintf("24. & High cost function intercepts  & %.2f & %.2f & %.2f & %.2f \\\\ \n", 
    ce_highera[1], ce_highera[2], ce_highera[3], ce_highera[4])
TabCE *= "\\midrule \n"
TabCE *= "\\end{tabular}"
write(dirtab * "TabCE.tex", TabCE)   



# ==============================================================
# Output single-number tex files for reference case
# ==============================================================
# Periods per year
open(dirsnt*"cal/T.tex", "w") do file
    write(file, string(p_ref.T))
end
# Months per period
open(dirsnt*"cal/monthsT.tex", "w") do file
    write(file, string(Int(round(12/p_ref.T, digits=0))))
end
# Number of years
open(dirsnt*"cal/maxY.tex", "w") do file
    write(file, string(Int(round(p_ref.maxY, digits=0))))
end

# Initial production rates, loop over regions
for (i, p) in enumerate([p_ref.q0[1], p_ref.q0[2], p_ref.q0[3], p_ref.q0[4]])
    open(dirsnt * "cal/q0_$(i).tex", "w") do file
        write(file, string(round(p, digits=1)))
    end
end
# Global Q0
open(dirsnt*"cal/Q0.tex", "w") do file
    write(file, string(round(sum(p_ref.q0), digits=1)))
end

# Initial price
open(dirsnt*"cal/P0.tex", "w") do file
    write(file, string(round(p_ref.Pref, digits=2)))
end
# Demand elasticity
open(dirsnt*"cal/delast.tex", "w") do file
    write(file, string(round(p_ref.delast, digits=1)))
end
# Ref case time for demand to go to zero
open(dirsnt*"cal/timetozero.tex", "w") do file
    write(file, string(Int(round(p_ref.timetozero, digits=0))))
end
# Initial demand growth rate
open(dirsnt*"cal/dgrowth.tex", "w") do file
    write(file, string(round(p_ref.dgr_ann*100, digits=2)))
end
open(dirsnt*"cal/dgrowthdec.tex", "w") do file
    write(file, string(round(1+p_ref.dgr_ann, digits=4)))
end

# Reserves, loop over regions
for (i, p) in enumerate([p_ref.x0[1], p_ref.x0[2], p_ref.x0[3], p_ref.x0[4]])
    open(dirsnt * "cal/x0_$(i).tex", "w") do file
        write(file, string(Int(round(p, digits=0))))
    end
end
# Global X0
open(dirsnt*"cal/X0.tex", "w") do file
    write(file, string(Int(round(sum(p_ref.x0), digits=0))))
end
# Reserves to production
open(dirsnt*"cal/RP0.tex", "w") do file
    write(file, string(Int(round(sum(p_ref.x0)/sum(p_ref.q0)*1000/365, digits=0))))
end

# Conventional decline rate
open(dirsnt*"cal/lambda1.tex", "w") do file
    write(file, string(Int(round(p_ref.lambday[1]*100, digits=0))))
end
# Shale decline rate
open(dirsnt*"cal/lambda4.tex", "w") do file
    write(file, string(Int(round(p_ref.lambday[4]*100, digits=0))))
end

# Core OPEC discount rate
open(dirsnt*"cal/r1.tex", "w") do file
    write(file, string(Int(round(p_ref.r_ann[1]*100, digits=0))))
end
# Others' discount rate
open(dirsnt*"cal/r4.tex", "w") do file
    write(file, string(Int(round(p_ref.r_ann[4]*100, digits=0))))
end
# Discount rate for emissions
open(dirsnt*"cal/rE.tex", "w") do file
    write(file, string(Int(round(p_ref.Er_ann*100, digits=0))))
end

# Cost intercepts in $/bbl, loop over regions
for i in 1:4
    open(dirsnt*"cal/alpha0_$(i).tex", "w") do file
        write(file, string(Int(round(p_ref.alpha[i], digits=0))))
    end
end
# Initial investment rates, loop over regions
for region in 1:4
    open(dirsnt * "cal/d0_$(region).tex", "w") do file
        write(file, string(round(d0[region], digits=2)))
    end
end
# Increase in $/bbl MC for a 1 mmbbl/d increase in capacity, loop over regions
for region in 1:4
    open(dirsnt * "cal/gamma_$(region).tex", "w") do file
        write(file, string(round(p_ref.gamma[region] / λ(m)[region], digits=2)))
    end
end
# Drilling cost elasticities, loop over regions
for region in 1:4
    open(dirsnt * "cal/ce_$(region).tex", "w") do file
        write(file, string(round(ce_ref[region], digits=2)))
    end
end
for region in 1:4
    write(dirsnt * "cal/ce_$(region).tex", @sprintf("%.2f", ce_ref[region]))
end
# Convergence parameters
write(dirsnt * "cal/gain.tex", @sprintf("%.1f", gain(m)))
write(dirsnt * "cal/ptol.tex", @sprintf("%.0e", pricetol(m)))
write(dirsnt * "cal/stol.tex", @sprintf("%.0e", stol(m)))
write(dirsnt * "cal/mutol.tex", @sprintf("%.0e", mutol(m)))


# ==============================================================
# Output single-number tex files for alternative specs
# ==============================================================
# Slow decline rate
open(dirsnt*"cal/lambda1slow.tex", "w") do file
    write(file, string(Int(round(p_dec06.lambday[1]*100, digits=0))))
end
# Global X0, Rystad 2PC
open(dirsnt*"cal/X0_lowres.tex", "w") do file
    write(file, string(Int(round(sum(p_lowres.x0), digits=0))))
end
# Global X0, IEA
open(dirsnt*"cal/X0_highres.tex", "w") do file
    write(file, string(Int(round(sum(p_highres.x0), digits=0))))
end
# Drilling cost elasticities for high cost intercept spec, loop over regions
for region in 1:4
    open(dirsnt * "cal/ce_highera_$(region).tex", "w") do file
        write(file, string(round(ce_highera[region], digits=2)))
    end
end


# ==============================================================
# Calibrated gammas, in units of $/bbl per mmbbl/d investment
# For copy/paste to README.md
# ==============================================================
lambda_ref = 1 .- (1 .- p_ref.lambday).^(1/p_ref.T)
gamma_ref = p_ref.gamma ./ lambda_ref
lambda_dec08 = 1 .- (1 .- p_dec08.lambday).^(1/p_dec08.T)
gamma_dec08 = p_dec08.gamma ./ lambda_dec08
lambda_dec06 = 1 .- (1 .- p_dec06.lambday).^(1/p_dec06.T)
gamma_dec06 = p_dec06.gamma ./ lambda_dec06
lambda_dec30 = 1 .- (1 .- p_dec30.lambday).^(1/p_dec30.T)
gamma_dec30 = p_dec30.gamma ./ lambda_dec30
gamma_n = p_n.gamma
gamma_highres = p_highres.gamma ./ lambda_ref
gamma_lowres = p_lowres.gamma ./ lambda_ref
gamma_dall09 = p_dall09.gamma ./ lambda_ref
gamma_dall03 = p_dall03.gamma ./ lambda_ref
gamma_n_lowres = p_n_lowres.gamma
gamma_n_dall03 = p_n_dall03.gamma
gamma_R = p_R.gamma ./ lambda_ref
gamma_nmp = p_nmp.gamma ./ lambda_ref
gamma_delhigh = p_delhigh.gamma ./ lambda_ref
gamma_dellow = p_dellow.gamma ./ lambda_ref
gamma_dng = p_dng.gamma ./ lambda_ref
gamma_dhg = p_dhg.gamma ./ lambda_ref
gamma_highera = p_highera.gamma ./ lambda_ref
