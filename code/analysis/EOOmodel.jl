# ==============================================================
# Model for End of Oil, non-AKS model
# ==============================================================


# ==============================================================
# Packages
# ==============================================================
# import Pkg; Pkg.add("Parameters")
# import Pkg; Pkg.add("Roots")
# import Pkg; Pkg.add("QuadGK")
# import Pkg; Pkg.add("NLsolve")
# import Pkg; Pkg.add("BenchmarkTools")
# import Pkg; Pkg.add("JLD2")
# import Pkg; Pkg.add("FileIO")
# import Pkg; Pkg.add("Printf")
# import Pkg; Pkg.add("Plots")

using Parameters, Roots, QuadGK, NLsolve, BenchmarkTools
using LinearAlgebra
using JLD2, FileIO
using Printf
using Plots


# ==============================================================
# Define abstract EOO type that will be the default for all methods
# Define EOO_N struct as a subtype of EOO
# ==============================================================
abstract type EOO end


# ==============================================================
# Define EOO_N struct as a subtype of EOO
# ==============================================================
@with_kw struct EOO_N <: EOO
    αin::Vector{Float64}        # MC at q = 0, in $/bbl
    γin::Vector{Float64}        # change in MC per change in Q ($/bbl per mmbbl/d)
    x0in::Vector{Float64}       # reserves in each region, billion bbl
    delast::Float64             # demand elasticity at Pref
    Pref::Float64               # reference oil price, $/bbl
    dref::Float64               # quantity demanded at Pref (mmbbl/d)
    timetozero::Float64         # years until demand for oil is zero
    shiftdelin::Float64 = 0.    # delay in years before demand decline starts
    drem::Float64 = 0.          # share of original demand remaining at end of demand decline
    dgr_ann::Float64 = 0.       # annual growth rate in demand
    dgrtin::Float64 = 7.        # number of years that demand grows
    r_ann::Vector{Float64}      # annual discount rate
    Er_ann::Float64             # annual discount rate for emissions discounting
    mp::Bool                    # market power flag
    maxY::Int64 = 200           # maximum years to simulate
    DY::Float64 = 365.25        # days per year
    T::Int64                    # periods per year (for discrete output and myopic model)
    stol::Float64 = 1e-4        # solver tolerance (on cumulative production)
    mutol::Float64 = 1e-4       # cutoff below which initial shadow values are considered zero
    maxit::Int64 = 100          # maximum number of iterations for solver
end

# Create accessor functions
fields = [:delast, :Pref, :dref, :timetozero, :drem, :dgr_ann,
     :r_ann, :Er_ann, :mp, :maxY, :DY, :T, :stol, :mutol, :maxit]
for field in fields
    @eval $field(m::EOO_N) = m.$field
end


# ==============================================================
# Short functions that define more model parameters
# ==============================================================
# Short functions for all EOO models
maxT(m::EOO) = maxY(m)*T(m)                         # Max number of periods to simulate
dslope(m::EOO) = delast(m) * dref(m) / Pref(m)      # demand slope, mmbbl/d per $/bbl
dint(m::EOO) = dref(m) - dslope(m) * Pref(m)        # demand intercept (Q in mmbbld/d when P=0)
dintP(m::EOO) = -dint(m) / dslope(m)                # inverse demand intercept
# Vector of discount factors for each period
Ediscvec(m::EOO) = [1.0 / ((1.0 + Er_ann(m)) ^ (t / T(m))) for t in 0:maxT(m)-1]

# Short functions specific to EOO_N
α(m::EOO_N) = m.αin                 # MC at q = 0, in $/bbl
γ(m::EOO_N) = m.γin                 # change in MC per change in Q ($/bbl per mmbbl/d)
x0(m::EOO_N) = m.x0in               # reserves in each region, billion bbl
shiftdel(m::EOO_N) = m.shiftdelin   # delay in years before demand decline starts
dgrt(m::EOO_N) = m.dgrtin           # number of years that demand grows
# annual shift down in demand, $/bbl
shift(m::EOO_N) = dintP(m) / (timetozero(m) - shiftdel(m))
# time at which demand shift stops  
shiftstop(m::EOO_N) = shiftdel(m) + (timetozero(m) - shiftdel(m)) * (1-drem(m))
r(m::EOO_N) = log.(1 .+ r_ann(m))       # Continuous time discount rate
dgr(m::EOO_N) = log.(1 .+ dgr_ann(m))   # Continuous time demand growth rate

# Number of regions (specific to EOO_N)
function NN(m::EOO_N)
    N = length(α(m))
    # Check that all region vectors have same lengths
    if length(γ(m)) != N || length(r_ann(m)) != N || length(x0(m)) != N
        error("Supply-side parameter input vectors are of different lengths")
    else
        return N
    end
end


# ==============================================================
# Marginal cost function (all EOO models)
# ==============================================================
function mc(m::EOO, q::Vector{Float64})
    # INPUTS
    # m: an instance of EOO
    # q: (Nx1) production in mmbbl/d

    # OUTPUTS
    # c: (Nx1) marginal cost in $/bbl

    c = α(m) .+ γ(m) .* q
    return c
end


# ==============================================================
# MC elasticity (all EOO models)
# ==============================================================
function mce(m::EOO, q::Vector{Float64})
    # INPUTS
    # m: an instance of EOO
    # q: (Nx1) production in mmbbl/d

    # OUTPUTS
    # ce: (Nx1) marginal cost elasticity

    c = mc(m, q)  # MC level
    ce = γ(m) .* q ./ c  # elasticity
    return ce
end   


# ==============================================================
# Quantity supplied as function of p and mu (EOO_N models)
# ==============================================================
function imc(m::EOO_N, p::Float64, mu::Vector{Float64})
    # INPUTS
    # m: an instance of EOO
    # p: price in $/bbl
    # mu: (Nx1) current shadow values, $/bbl

    # OUTPUTS
    # q: (Nx1) production in mmbbl/d

    num = max.(0, p .- α(m) .- mu)  # numerator
    q = num ./ γ(m)

    # FOC for region 1 if it exerts market power
    if mp(m)
        q[1] = num[1] / (γ(m)[1] - 1/dslope(m))
    end
    return q
end


# ==============================================================
# Quantity demanded (all EOO models)
# ==============================================================
function dem(m::EOO, p::Float64, t::Vector{Float64}, 
    td0::Real, decreaseflag::Bool)
    # INPUTS
    # m: an instance of EOO
    # p: price in $/bbl
    # t: time; shifts demand curve (vector when called from an AKS model)
    # td0: time that sets level of demand decline for all t if decreaseflag==false
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # Note that time units depend on the instance of EOO. 
    # EOO_N uses years. EOO_AKS uses periods

    # OUTPUTS
    # q: quantity demanded in mmbbl/d

    # Compute time to multiply against shift(m)
    if decreaseflag==true
        ts = min.(max.(0, t .- shiftdel(m)), shiftstop(m) - shiftdel(m))
    else
        ts = min(max(0, td0 - shiftdel(m)), shiftstop(m) - shiftdel(m))
    end
    
    # Quantity demanded ignoring demand growth
    qng = max.(0, dint(m) .+ dslope(m) * (p .+ shift(m) * ts))

    # Account for demand growth
    q = qng .* exp.(dgr(m) .* min.(t, dgrt(m)))
    return q
end


# ==============================================================
# Inverse demand (all EOO models)
# ==============================================================
function idem(m::EOO, q::Vector{Float64}, t::Vector{Float64}, 
    td0::Real, decreaseflag::Bool)
    # INPUTS
    # m: an instance of EOO
    # q: quantity demanded in mmbbl/d (vector when called from an AKS model)
    # t: time in years; shifts demand curve (vector when called from an AKS model)
    # td0: time that sets level of demand decline for all t if decreaseflag==false
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # Note that time units depend on the instance of EOO. 
    # EOO_N uses years. EOO_AKS uses periods

    # OUTPUTS
    # p: price in $/bbl

    # First reduce q to degrowthed quantity
    qng = q ./ exp.(dgr(m) .* min.(t, dgrt(m)))

    # Compute time to multiply against shift(m)
    if decreaseflag==true
        ts = min.(max.(0, t .- shiftdel(m)), shiftstop(m) - shiftdel(m))
    else
        ts = min(max(0, td0 - shiftdel(m)), shiftstop(m) - shiftdel(m))
    end

    p = max.(0, (qng .- dint(m)) / dslope(m) .- shift(m) * ts)
    p = min.(dintP(m), p)  # limit to inverse demand intercept
    return p
end


# ==============================================================
# Diff between quantity demanded and supplied at a given p and mus (EOO_N models)
# ==============================================================
function qdiff(m::EOO_N, p::Float64, mu::Vector{Float64}, 
    t::Float64, td0::Float64, decreaseflag::Bool)
    # INPUTS
    # m: an instance of EOO_N
    # p: candidate price, $/bbl
    # mu: (Nx1) current shadow values, $/bbl
    # t: time in years, sets demand
    # td0: Sets level of demand decline if decreaseflag==false
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases

    # OUTPUTS
    # d: (Nx1) diff between quantities demanded and supplied

    qd = dem(m, p, [t], td0, decreaseflag)      # quantity demanded
    qs = imc(m, p, mu)     # quantities supplied
    d = qd .- sum(qs)       # diff
    return d
end
    

# ==============================================================
# Solve for q and p given current shadow values and time (EOO_N models)
# ==============================================================
function solvepq(m::EOO_N, mu::Vector{Float64}, 
    t::Float64, td0::Float64, decreaseflag::Bool)
    # INPUTS
    # m: an instance of EOO_N
    # mu: (Nx1) current shadow values, $/bbl
    # t: time in years
    # td0: Sets level of demand decline if decreaseflag==false
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases

    # OUTPUTS
    # q: (Nx1) quantities in mmbbl/d
    # p: price, $/bbl

    pguess = minimum(α(m) .+ mu)        # Guess for p
    # Solve for p given guess pguess
    p = fzero(x -> qdiff(m, x, mu, t, td0, decreaseflag), pguess, xtol=1e-14)
    p = min(dintP(m), p)  # limit to inverse demand intercept
    # Compute quantities
    q = imc(m, p, mu)

    return q, p
end


# ==============================================================
# Compute cumulative production through T given initial shadow value (EOO_N models)
# ==============================================================
function cumprod(m::EOO_N, Te::Float64, mu0::Vector{Float64}, 
    decreaseflag::Bool, td0::Real=0.0)
    # INPUTS
    # m: an instance of EOO_N
    # Te: time in years at which to stop cumulating (possibly infinite)
    # mu0 (Nx1): initial shadow values (at t=0), $/bbl
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: (optional): initial time in years (positions the initial demand curve)

    # OUTPUTS
    # Q: (Nx1) cumulative production in billion bbl

    # max number of years to cumulate
    tmax = min(Te, maxY(m))

    # Check that choke price at td0 exceeds mc + mu0 for at least one region
    # Account for maximum possible future demand growth
    chokeP = idem(m, [0.], [td0], td0, decreaseflag) * exp(dgr(m) * dgrt(m))
    if maximum(chokeP .- α(m) .- mu0 .* exp.(r(m) .* td0)) <= 0 || Te <= td0
        Q = zeros(Float64, NN(m))              # no production
    else                    
        # First define utility function to assist integration
        function util(m::EOO_N, t::Float64, mu0::Vector{Float64}, td0::Float64, 
            decreaseflag::Bool)
            q, _ = solvepq(m, mu0 .* exp.(r(m) .* t), t, td0, decreaseflag)
            return q
        end
        Q, _ = quadgk(t -> util(m, t, mu0, td0, decreaseflag), td0, tmax, rtol=1e-12)
        Q *= DY(m) / 1000        # convert to billions of bbl / year
    end
    return Q
end


# ==============================================================
# Solve for optimal production (EOO_N models)
# ==============================================================
function optpath(m::EOO_N, decreaseflag::Bool, td0::Real, x0::Vector{Float64}, 
    mu0guess::Vector{Float64}=fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO_N
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: initial time in years (positions the initial demand curve)
    # x0: (Nx1) initial reserves, in billions of bbl
    # mu0guess (optional): (Nx1) guess of initial shadow value in $/bbl

    # OUTPUTS
    # Q: (Nx1) cumulative production in billion bbl
    # pvec: (Tx1) annual time series of prices in $/bbl
    # qvec: (TxN) annual time series of production rates in mmbbl/d
    # mu0: (Nx1) initial (t=0) shadow value in $/bbl
    # muvec: (TxN) time series of current shadow values in $/bbl
    # Te (Nx1) approx. times in years when production stops for each region

    # Note: for all indices of mu0guess that are zero, the reported mu0 will be zero

    # Check if Q < reserves at mu0 = 0 for all fields
    Q = cumprod(m, Inf, zeros(NN(m)), decreaseflag, td0)
    if all(Q .<= x0)
        mu0 = zeros(NN(m))
    elseif all(x -> x == 0. || isinf(x), mu0guess)
        # mu0guess is all zero or Inf; don't search
        mu0 = mu0guess
    elseif x0 == zeros(NN(m))
        Q = zeros(NN(m))
        pvec = zeros(maxT(m))
        qvec = zeros(maxT(m), NN(m))
        mu0 = zeros(NN(m))
    else
        # Need to search for mu0
        # Define utility function for the complementarity problem, lx is log(mu0)
        function utilc(m::EOO_N, lx::Vector{Float64}, decreaseflag::Bool, 
            td0::Real, x0::Vector{Float64})
            cumD = cumprod(m, Inf, exp.(lx), decreaseflag, td0) 
            rawD = x0 - cumD    # raw diff between reserves and drilling            
            # adjustment -- smooth toward complementary slackness at mu=0
            adj = ones(NN(m))
            adjind = (rawD .> 0.) .& (exp.(lx).<mutol(m))
            adj[adjind] .= exp.(lx[adjind]).^3 / mutol(m)^3
            DD = rawD .* adj
            # second adjustment -- penalize high mu's with zero drilling
            # helps prevent algorithm from getting stuck on high mu
            adj2 = zeros(NN(m))
            adj2[cumD .== 0] .= exp.(lx[cumD .== 0]) * stol(m) / mutol(m) * 10
            DD2 = DD + adj2.^2
            return DD2
        end
        # Final outer utility function that enables iteration over only non-zero / non-inf mu's
        nzind = mu0guess .> 0 .&& mu0guess .< Inf
        infind = mu0guess .== Inf
        function utilfm(m::EOO_N, lx::Vector{Float64}, decreaseflag::Bool, 
            td0::Real, x0::Vector{Float64}, nzind::BitVector, infind::BitVector)
            lxall = ones(NN(m))*(-20)  # initialize with small values
            lxall[infind] .= Inf
            lxall[nzind] = lx
            Dall = utilc(m, lxall, decreaseflag, td0, x0)
            return Dall[nzind]
        end
        # Solve for mu0
        result = nlsolve(lx -> utilfm(m, lx, decreaseflag, td0, x0, nzind, infind), 
            log.(mu0guess[nzind]), ftol=stol(m), show_trace=true, iterations=maxit(m))
        mu0 = zeros(NN(m))
        mu0[mu0guess.==Inf] .= Inf
        mu0[nzind] = exp.(result.zero)
        # Replace small values with zero
        mu0[mu0 .< mutol(m)] .= 0.
        # Get extraction and prices at mu0
        Q = cumprod(m, Inf, mu0, decreaseflag, td0)
    end
    # Create time series outputs. Time measured in periods, not years
    pvec = zeros(maxT(m))
    qvec = zeros(maxT(m), NN(m))
    muvec = zeros(maxT(m), NN(m))
    td0Tint = Int(floor(td0*T(m))) + 1   # handle non-integer td0*T
    for t in td0Tint:maxT(m)
        mu = mu0 .* exp.(r(m) .* (td0 + (t-td0Tint) / T(m)))
        q, p = solvepq(m, mu, (t-1) / T(m), td0, decreaseflag)
        muvec[t, :] = mu'
        pvec[t] = p
        qvec[t, :] = q'
    end
    # Stop times in years
    Te = zeros(NN(m))
    for n in 1:NN(m)
        if Q[n] == 0.
            Te[n] = 0.
        else
            nz = findall(qvec[:, n] .> 0)
            Te[n] = maximum(nz) / T(m)
        end
    end
    return Q, pvec, qvec, mu0, muvec, Te
end
   

# ==============================================================
# Myopic path in which firms don't believe demand is declining (EOO_N models)
# ==============================================================
function myopicpath(m::EOO_N, mu0guess::Vector{Float64} = fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO_N
    # mu0guess: (optional) (Nx1) guess of initial shadow value in $ per mmbbl/d

    # OUTPUTS
    # Q: (Nx1) cumulative production in billion bbl
    # pvec: (Tx1) time series of prices in $/bbl
    # qvec (TxN): time series of production rates in mmbbl/d
    # muvec: (TxN) time series of initial shadow values in $/bbl
    
    # Initialize loop
    t = 0
    endflag = false
    pvec = zeros(maxT(m))
    qvec = zeros(maxT(m), NN(m))
    Q = zeros(NN(m))
    muvec = zeros(maxT(m), NN(m))
    xt = x0(m)  # initial reserves
    lastp = min(shiftstop(m)*T(m), maxT(m)-1)   # last period of loop

    # Loop over periods
    while t <= lastp && endflag == false
        # Solve model assuming demand is fixed at its level halfway through the period
        _, pt, qt, mut, _, _ = optpath(m, false, (t+0.5) / T(m), xt, mu0guess)
        # Compute one period of cumulative production on this path
        Qt = cumprod(m, (t + 1.5) / T(m), mut, false, (t+0.5) / T(m))
        # Increment time and store results
        t += 1
        pvec[t] = pt[t]
        qvec[t, :] = qt[t, :]
        muvec[t, :] = mut'
        # Decrease reserves
        xt -= Qt
        println("Current t and value of xt: ", [t-1 xt'])
        mu0guess = mut  # update guess
        mu0guess[xt .<= 0.] .= Inf  # set mu0guess to Inf if reserves are exhausted
        # Improve guess
        guessratio = ones(NN(m))
        if t>=2
            for i = 1:NN(m)
                if muvec[t-1, i] > 0. && muvec[t, i] < Inf
                    guessratio[i] = mut[i] / muvec[t-1, i]
                end
            end
        end
        mu0guess = mu0guess .* guessratio
        println("Current mut and mu0guess: ", [mut' mu0guess'])
        # Trigger endflag if reserves are exhausted
        if maximum(xt) <= 0.
            endflag = true
        end
    end

    # Final set of periods when demand is at steady state again (possibly zero)
    if lastp==maxT(m)-1 || drem(m)==0. || maximum(xt) <= 0.  
        Qe = zeros(NN(m))  # nothing else to do; all remaining periods have zero q and p
    else
        Qe, pte, qte, mute, _, _ = optpath(m, false, (t+0.5) / T(m), xt, mu0guess)
        pvec[t+1:end] = pte[t+1:end]
        qvec[t+1:end, :] = qte[t+1:end, :]
        muvec[t+1:end, :] = repeat(mute', outer=[length(t+1:maxT(m)), 1])
    end

    # Cumulative production
    Q = x0(m) - xt + Qe
    return Q, pvec, qvec, muvec
end


# ==============================================================
# Solve for PV and profits from solved model (EOO_N models)
# ==============================================================
function cumprofits(m::EOO_N, pvec::Vector{Float64}, qvec::Matrix{Float64})
    # INPUTS
    # m: an instance of EOO_N
    # pvec: (Tx1) annual time series of prices in $/bbl
    # qvec (TxN): rate of production in each period, mmbbl/d

    # OUTPUTS
    # PVpi: (Nx1) PDV of profits for each region ($billion)
    # pivec: (TxN) time series of instantaneous profit flows ($billion / period)

    # Compute per-period profits (TxN)
    pivec = (pvec .* qvec - α(m)' .* qvec - γ(m)' .* qvec.^2 / 2) * DY(m) / 1000 / T(m) 

    # Compute PV
    Tvec = 0:maxT(m)-1
    Deltavec = (1 ./ (1 .+ r_ann(m)))'.^(Tvec/T(m))   # discount factor for each period
    Discpivec = pivec .* Deltavec
    PVpi = vec(sum(Discpivec, dims=1))

    return PVpi, pivec
end


# ==============================================================
# Compute first period supply elasticity (EOO_N models)
# ==============================================================
function supplyelast(m::EOO_N, mu0guess::Vector{Float64}=fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO_N

    # OUTPUTS
    # elasttot: total supply elasticity
    # elast: supply elasticity for each region
    
    # Outcomes at actual first period demand
    _, pvec, qvec, _, _, _ = optpath(m, false, 0., x0(m), mu0guess)

    # Outcomes at shifted first year demand
    _, pvecs, qvecs, _, _, _ = optpath(m, false, 1. / T(m), x0(m), mu0guess)

    elast = ((qvecs[2, :] .- qvec[1, :]) ./ (pvecs[2] .- pvec[1]) 
        ./ (qvecs[2, :] .+ qvec[1, :]) .* (pvecs[2] .+ pvec[1]))

    elasttot = ((sum(qvecs[2, :]) - sum(qvec[1, :])) / (pvecs[2] - pvec[1]) 
        / (sum(qvecs[2, :]) + sum(qvec[1, :])) * (pvecs[2] + pvec[1]))
    return elasttot, elast
end


# ==============================================================
# Summarize quantities and profits (EOO_N models)
# ==============================================================
function sum_Q_pi(m::EOO_N, Q, qvec::Matrix{Float64}, pvec::Vector{Float64})
    # INPUTS
    # m: an instance of EOO
    # Q: (Nx1) cumulative production in billion bbl
    # qvec: (TxN) time series of production rates in mmbbl/d
    # pvec: (Tx1) annual time series of prices in $/bbl

    # OUTPUTS
    # elasttot: total supply elasticity
    # elast: supply elasticity for each region
    
    # Summaries of production
    QT = sum(Q)         # total production in billion bbl
    # Discounted production by region, billion bbl
    PVQ = qvec' * Ediscvec(m) * DY(m) / 1000 / T(m)
    PVQT = sum(PVQ)     # total discounted production in billion bbl

    # PV of profits and per period profits, by region
    PVpi, pivec = cumprofits(m, pvec, qvec)
    
    return QT, PVQT, PVQ, PVpi, pivec
end


# ==============================================================
# Return diff between q0 and first period production
# ==============================================================
function q0diff(m::EOO_N, q0in::Vector{Float64}, γguess::Vector{Float64},
    mu0guess::Vector{Float64}=fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO_N
    # q0in: (Nx1) actual first period production in mmbbl/d
    # γguess: (Nx1) guess of γ parameter vector, in $/bbl per mmbbl/d
    # mu0guess (optional): (Nx1) guess of initial shadow value in $/bbl

    # OUTPUTS
    # diff: q0in minus first period production

    # Instantiate model with guessed γ
    mγ = EOO_N(αin=m.αin, γin=γguess, x0in=m.x0in, delast=m.delast, Pref=m.Pref, 
        dref=m.dref, timetozero=m.timetozero, shiftdelin=m.shiftdelin, drem=m.drem, 
        dgr_ann=m.dgr_ann, r_ann=m.r_ann, Er_ann=m.Er_ann, mp=m.mp, maxY=m.maxY, 
        DY=m.DY, T=m.T, stol=m.stol, mutol=m.mutol, maxit=m.maxit)
    
    # Simulate model and obtain production
    _, _, qvec, _, _, _ = optpath(mγ, false, 0., x0(mγ), mu0guess)

    # Outcomes at shifted first year demand
    diff = q0in - vec(qvec[1, :])
    return diff
end


# ==============================================================
# Solve for γ that equates q0 and first period production
# ==============================================================
function findgamma(m::EOO, q0in::Vector{Float64}, γguess::Vector{Float64},
    mu0guess::Vector{Float64}=fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO
    # q0in: (Nx1) actual first period production in mmbbl/d
    # γguess: (Nx1) guess of γ parameter vector
    # mu0guess (optional): (Nx1) guess of initial shadow value

    # OUTPUTS
    # γout: (Nx1) γ parameter vector that equates q0in and first period production

    # Solve for γ
    result = nlsolve(x -> q0diff(m, q0in, x, mu0guess), 
            γguess, ftol=stol(m), show_trace=true, iterations=maxit(m))
    γout = result.zero
    return γout
end


# ==============================================================
# Function to run non-AKS models given input parameter tuple
# ==============================================================
function runnoonAKSmodel(pin::NamedTuple, unant::Bool = true)
    # INPUTS
    # pin: tuple of input parameters
    # unant: flag for whether to run unanticipated demand decline

    # OUTPUTS
    # Oi: tuple of baseline outputs
    # O: tuple of anticipated demand decline outputs
    # Om: tuple of myopic demand decline outputs
    # Op: tuple of percent changes outputs
    m = EOO_N(αin=pin.alpha, γin=pin.gamma, x0in=pin.x0, 
        delast=pin.delast, Pref=pin.Pref, dref=pin.dref, 
        dgr_ann=pin.dgr_ann, timetozero=pin.timetozero, 
        r_ann=pin.r_ann, Er_ann=pin.Er_ann, mp=pin.mp, T=pin.T,
        shiftdelin=pin.shiftdelin, drem=pin.drem, maxY=pin.maxY, mutol=pin.mutol)
    Qi, pveci, qveci, mu0i, muveci, Tei = optpath(m, false, 0., x0(m))
    mci_bbl = mc(m, vec(qveci[1,:]))             # initial MC in $/bbl
    QTi, PVQTi, PVQi, PVpii, piveci = sum_Q_pi(m, Qi, qveci, pveci)
    # Demand decline   
    Qd, pvecd, qvecd, mu0d, muvecd, Ted = optpath(m, true, 0., x0(m))
    QTd, PVQTd, PVQd, PVpid, pivecd = sum_Q_pi(m, Qd, qvecd, pvecd)
    # Myopic
    if unant==true
        Qm, pvecm, qvecm, muvecm = myopicpath(m, mu0i)
        QTm, PVQTm, PVQm, PVpim, pivecm = sum_Q_pi(m, Qm, qvecm, pvecm)
    else
        Qm = ones(NN(m)); pvecm = ones(maxT(m)); qvecm = ones(maxT(m),NN(m)); muvecm = ones(maxT(m),NN(m));
        QTm = 1; PVQTm = 1; PVQm = ones(NN(m)); 
        PVpim = ones(NN(m)); pivecm = ones(maxT(m),NN(m));
    end
    # Percent changes
    PctIncQdQiTot = (QTd-QTi) / QTi * 100; PctIncPVQdQiTot = (PVQTd-PVQTi) / PVQTi * 100
    if unant==true
        PctIncQdQmTot = (QTd-QTm) / QTm * 100; PctIncPVQdQmTot = (PVQTd-PVQTm) / PVQTm * 100
        PctIncQmQiTot = (QTm-QTi) / QTi * 100; PctIncPVQmQiTot = (PVQTm-PVQTi) / PVQTi * 100
    else
        PctIncQdQmTot = 1; PctIncPVQdQmTot = 1; PctIncQmQiTot = 1; PctIncPVQmQiTot = 1
    end
    # Create output touples
    Oi = (Q = Qi, pvec = pveci, qvec = qveci, 
        mu0 = mu0i, muvec = muveci, Te = Tei,
        mc_bbl = mci_bbl, QT = QTi, 
        PVQT = PVQTi, PVQ = PVQi, PVpi = PVpii, pivec = piveci)
    Od = (Q = Qd, pvec = pvecd, qvec = qvecd,
        mu0 = mu0d, muvec = muvecd, Te = Ted, QT = QTd, 
        PVQT = PVQTd, PVQ = PVQd, PVpi = PVpid, pivec = pivecd)
    Om = (Q = Qm, pvec = pvecm, qvec = qvecm,
        muvec = muvecm, QT = QTm, PVQT = PVQTm, PVQ = PVQm, 
        PVpi = PVpim, pivec = pivecm)
    Op = (PctIncQdQmTot = PctIncQdQmTot,
        PctIncPVQdQmTot = PctIncPVQdQmTot, PctIncQdQiTot = PctIncQdQiTot,
        PctIncPVQdQiTot = PctIncPVQdQiTot, PctIncQmQiTot = PctIncQmQiTot,
        PctIncPVQmQiTot = PctIncPVQmQiTot)
    return Oi, Od, Om, Op
end


# ==============================================================
# Function for users to calibrate and run entire model
# ==============================================================
function runall(pin, p_AKS, p_estslopes, p_unant, mu0_bbl_guess)
    # INPUTS
    # pin: tuple of input parameters
    # p_AKS: flag for whether to run AKS investment model
    # p_estslopes: flag for whether to calibrate cost function slopes
    # p_unant: flag for whether to run unanticipated demand decline
    # mu0_bbl_guess: (Nx1) guess of initial shadow value in $/bbl

    # OUTPUTS

    # Calibrate cost function slopes if requested
    if p_estslopes==true
        # Instantiate non-AKS model with guessed γ
        m_n_guess = EOO_N(αin=pin.alpha, γin=pin.gamma, x0in=pin.x0, 
            delast=pin.delast, Pref=pin.Pref, dref=pin.dref, 
            dgr_ann=pin.dgr_ann, timetozero=pin.timetozero, 
            r_ann=pin.r_ann, Er_ann=pin.Er_ann, mp=pin.mp, T=pin.T,
            shiftdelin=pin.shiftdelin, drem=pin.drem, maxY=pin.maxY, mutol=pin.mutol)
        if p_AKS==true
            # Run non-AKS model to get price vector guess
            _, pveci_n, ~, mu0i_n, _, _ = optpath(m_n_guess, false, 0., x0(m_n_guess))
            # Instantiate AKS model and solve for gamma
            m_guess = EOO_AKS(eoo_n=m_n_guess, lambday=pin.lambday, q0=pin.q0, pguess=pveci_n)
            p_gamma = findgamma(m_guess, pin.q0, pin.gamma, mu0_bbl_guess .* RRC(m_guess))
            # Store in touple
            p = merge(pin, Dict(:gamma => p_gamma))
        else
            # Solve for gamma
            p_gamma_n = findgamma(m_n_guess, pin.q0, pin.gamma, mu0_bbl_guess)
            # Store in touple
            p = merge(pin, Dict(:gamma => p_gamma_n))
        end
    else
        p = pin
    end

    # Run models
    if p_AKS==true
        Baseline, Ant, Unant, PctInc = runAKSmodel(p, p_unant)
    else
        Baseline, Ant, Unant, PctInc = runnoonAKSmodel(p, p_unant)
    end

    return Baseline, Ant, Unant, PctInc
end








