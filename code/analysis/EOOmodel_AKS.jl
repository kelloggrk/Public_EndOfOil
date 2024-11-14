# ==============================================================
# Model for End of Oil, with AKS-style investment
# ==============================================================

# ==============================================================
# Define EOO_AKS struct as a subtype of abstract type EOO
# ==============================================================

@with_kw struct EOO_AKS <: EOO
    eoo_n::EOO_N                # EOO_N model struct
    lambday::Vector{Float64}    # (Nx1) annual decline rates
    q0::Vector{Float64}         # (Nx1) initial production rates mmbbl/d
    pguess::Vector{Float64}     # (maxTx1) guess of time path of prices
    gain::Float64 = 0.2         # gain when iterating on price path 
    maxPit::Int64 = 200         # maximum number of iterations on price path search
    pricetol::Float64 = 1e-5    # tolerance for price path search
end

# Forward accessor functions from EOO_N
for method in (:delast, :Pref, :dref, :timetozero, :drem, :dgr_ann,
    :r_ann, :Er_ann, :mp, :maxY, :DY, :T, :stol, :mutol, :maxit)
    @eval $method(m::EOO_AKS) = $method(m.eoo_n)
end

# New accessor functions
lambday(m::EOO_AKS) = m.lambday
q0(m::EOO_AKS) = m.q0
pguess(m::EOO_AKS) = m.pguess
gain(m::EOO_AKS) = m.gain
maxPit(m::EOO_AKS) = m.maxPit
pricetol(m::EOO_AKS) = m.pricetol


# ==============================================================
# Short functions that define more model parameters
# ==============================================================
r(m::EOO_AKS) = (1 .+ r_ann(m)).^(1/T(m)) .- 1      # per-period discount rate
δ(m::EOO_AKS) = 1 ./ (1 .+ r(m))                    # per-period discount factor
λ(m::EOO_AKS) = 1 .- (1 .- lambday(m)).^(1/T(m))    # per-period decline rate
dgr_p(m::EOO_AKS) = (1 .+ dgr_ann(m)).^(1/T(m)) .- 1    # per-period demand growth rate
dgr(m::EOO_AKS) = log.(1 .+ dgr_p(m))   # Continuous time demand growth rate, in period^-1
dgrt(m::EOO_AKS) = dgrt(m.eoo_n) * T(m)             # number of periods of demand growth

# Reserves (billion bbl) per mmbbl/d of capacity
RC(m::EOO_AKS) = DY(m) / 1000 / T(m) ./ λ(m);
# Time discounted reserves (billion bbl) per mmbbl/d of capacity
# Accounts for fact that decline starts two periods after drilling
RRC(m::EOO_AKS) = DY(m) / 1000 / T(m) ./ (r(m) .+ λ(m));

# Convert costs from oilmodelN to capacity investment costs ($billion per mmbbl/d invested)
α(m::EOO_AKS) = α(m.eoo_n) .* RRC(m);              # intercept of MC
γ(m::EOO_AKS) = γ(m.eoo_n) .* RRC(m) ./ λ(m);      # slope ($b per mmbbl/d per mmbbl/d invested)

# delay before demand decline starts, in periods
shiftdel(m::EOO_AKS) = shiftdel(m.eoo_n) * T(m)
shift(m::EOO_AKS) = shift(m.eoo_n) / T(m)   # per-period shift down in demand
# time at which demand shift stops  
shiftstop(m::EOO_AKS) = shiftdel(m) + (timetozero(m)*T(m) - shiftdel(m)) * (1-drem(m))

# Convert reserves to remaining mmbbl/d of capacity to drill
# Subtract off reserves produced through initial capacity q0
x0(m::EOO_AKS) = max.(0., x0(m.eoo_n) ./ RC(m) .- q0(m));

# Number of regions (specific to EOO_AKS)
function NN(m::EOO_AKS)
    N = length(α(m))
    # Check that all region vectors have same lengths
    if (length(γ(m)) != N || length(r_ann(m)) != N || length(x0(m)) != N 
        || length(lambday(m)) != N || length(q0(m)) != N)
        error("Supply-side parameter input vectors are of different lengths")
    else
        return N
    end
end


# ==============================================================
# Quantity supplied as function of theta and mu
# ==============================================================
function imc(m::EOO_AKS, θ::Vector{Float64}, mu::Vector{Float64})
    # INPUTS
    # m: an instance of EOO_AKS
    # θ: (Nx1) thetas in $/bbl
    # mu: (Nx1) current shadow values, $/bbl

    # OUTPUTS
    # d: (Nx1) rate of drilling (capacity investment), mmbbl/d per period

    # Note: market power accounted for in ctheta for EOO_AKS models

    num = max.(0, θ .- α(m) .- mu)  # numerator
    d = num ./ γ(m)
    return d
end


# ==============================================================
# Compute theta given a price vector
# ==============================================================
function ctheta(m::EOO_AKS, pvec::Vector{Float64}, q1vec::Vector{Float64})
    # INPUTS
    # m: an instance of EOO_AKS
    # pvec: vector of prices in $/bbl. Starts after the period corresponding to theta
    # q1vec: vector of region 1's quantities

    # OUTPUTS
    # theta: value of capacity in $billion per mmbbl/d

    Tp = length(pvec)
    Tvec = 1:Tp
    d = (1 .- λ(m)) .* δ(m)     # discount factors accounting for both decline and interest
    dvec = d .^ Tvec'           # vector of discount factors

    # Create marginal revenue matrix. This is just price for all but region 1
    mrmat = repeat(pvec, 1, NN(m))            # Repeat pvec for N columns
    # Change MR for region 1 if it has market power
    if mp(m)
        mrmat[:, 1] = mrmat[:, 1] .+ 1 / dslope(m) * q1vec
        mrmat[:, 1] = max.(0, mrmat[:, 1])
    end
    
    # Build theta by looping over regions
    θ = zeros(NN(m))
    for i in 1:NN(m)
        # Sum to get theta over length of mrmat
        # Divide by 1-lambda since decline doesn't start until 2nd period, but time
        # discounting starts in period 1
        θ[i] = DY(m) / 1000 / T(m) / (1 - λ(m)[i]) * dot(dvec[i, :], mrmat[:, i])
        # Extrapolate beyond pvec, assuming constant price after end of vector
        θe = mrmat[Tp, i] * d[i]^(Tp + 1) / (1 - d[i]) * DY(m) / 1000 / T(m) / (1 - λ(m)[i])
        θ[i] += θe
    end
    return θ
end


# ==============================================================
# Compute time series of drilling given a price vector and initial shadow vals
# ==============================================================
function drillpath(m::EOO_AKS, pvec::Vector{Float64}, 
    q1vec::Vector{Float64}, mu::Vector{Float64})
    # INPUTS
    # m: an instance of EOO_AKS
    # pvec: (Tpx1) vector of prices in $/bbl
    # mu: (Nx1) initial shadow values, $/bbl
    # q1vec: (Tpx1) vector of region 1's quantities

    # OUTPUTS
    # dvec (TpxN): rate of drilling (capacity investment) in each period, mmbbl/d per period

    # Initialize output
    Tp = length(pvec)
    dvec = zeros(Tp, NN(m))

    # Loop over time periods
    for t = 1:(Tp - 1)
        θ = ctheta(m, pvec[(t + 1):Tp], q1vec[(t + 1):Tp])  # compute theta at t
        mut = mu .* (1 .+ r(m)) .^ (t - 1)                  # current shadow val
        at = imc(m, θ, mut)                                 # current drilling
        dvec[t, :] = at'
    end
    return dvec
end


# ==============================================================
# Compute price and quantity paths given drilling paths and initial production
# ==============================================================
function pqpath(m::EOO_AKS, dvec::Matrix{Float64}, q0::Vector{Float64}, 
    decreaseflag::Bool, td0::Int64)
    # INPUTS
    # m: an instance of EOO_AKS
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period
    # q0: (Nx1) production in first period
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: initial time period (=1 for period when prod = q0)

    # OUTPUTS
    # pvec (Tx1): price in each period, $/bbl
    # qvec (TxN): rate of production in each period, mmbbl/d

    # Compute quantity paths
    qvec = zeros(maxT(m), NN(m))
    qvec[1, :] = q0'
    for t = 2:maxT(m)
        qvec[t, :] = qvec[t - 1, :]' .* (1 .- λ(m)') .+ dvec[t - 1, :]'
    end

    # Compute price paths
    Q = vec(sum(qvec, dims=2))      # total production each period
    Tvec = 0.:(maxT(m) - 1.)        # times for demand shift
    times = Float64[td0 - 1 + t for t in Tvec]
    pvec = idem(m, Q, times, td0-1, decreaseflag)        # eqbm prices

    # If price is zero, pro-rate down quantities so that total Q does not exceed
    # quantity demanded
    dint = dem(m, 0., times, td0-1, decreaseflag)
    p0ind = findall(x -> x == 0, pvec)
    if !isempty(p0ind)
        for t in p0ind
            Qt = sum(qvec[t, :])
            if Qt > 0
                ratio = dint[t] / Qt
                qvec[t, :] .*= ratio            
            end
        end
    end
    return pvec, qvec
end


# ==============================================================
# Compute time series of prices, drilling, and production given 
# initial shadow vals and initial production
# ==============================================================
function solvepq(m::EOO_AKS, mu0::Vector{Float64}, q0::Vector{Float64}, 
    decreaseflag::Bool, td0::Int64, pvecg::Vector{Float64} = pguess(m))
    # INPUTS
    # m: an instance of EOO_AKS
    # mu0: (Nx1) initial (first period) shadow values, $/bbl
    # q0: (Nx1) production in first period
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: initial time period (positions the initial demand curve)
    # pvecg: (optional) initial price vector guess

    # OUTPUTS
    # pvec (Tx1): price in each period, $/bbl
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period
    # qvec (TxN): rate of production in each period, mmbbl/d

    # initial guess of q1 each period is declined initial production
    Tvec = 0:(maxT(m) - 1)
    q1vecg = q0[1] * (1 .- λ(m)[1]) .^ Tvec

    # Iterate to convergence
    Change = Inf; pvec = zeros(NN(m)); 
    qvec = zeros(maxT(m), NN(m)); dvec = zeros(maxT(m), NN(m))
    for i = 1:maxPit(m)
        dvec = drillpath(m, pvecg, q1vecg, mu0)      # get drilling path
        # Compute price path consistent with dvec
        pvec, qvec = pqpath(m, dvec, q0, decreaseflag, td0)
        # Update guess
        pvecg = pvecg .+ gain(m) * (pvec .- pvecg)
        q1vecg = q1vecg .+ gain(m) * (qvec[:, 1] .- q1vecg)
        # Check convergence
        Change = maximum(abs.(pvec .- pvecg))
        if Change <= pricetol(m)
            break
        end
    end
    # Tell user what happened
    if Change > pricetol(m)
        @warn "Failure to converge on eqbm price series"
    end
    return pvec, dvec, qvec
end


# ==============================================================
# Compute cumulative drilling (through all time),
# given initial prod and initial shadow value
# ==============================================================
function cumdrill(m::EOO_AKS, mu0::Vector{Float64}, q0::Vector{Float64}, 
    decreaseflag::Bool, td0::Int64, pvecg::Vector{Float64} = pguess(m))
    # INPUTS
    # m: an instance of EOO_AKS
    # mu0: (Nx1) initial (first period) shadow values, $/bbl
    # q0: (Nx1) production in first period
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: initial time period (positions the initial demand curve)
    # pvecg: (optional) initial price vector guess

    # OUTPUTS
    # D: (Nx1) cumulative drilling capacity added in mmbbl/d

    # Obtain vectors of prices, drilling, and production
    _, dvec, _ = solvepq(m, mu0, q0, decreaseflag, td0, pvecg)

    # Cumulate drilling
    D = vec(sum(dvec, dims=1))
    return D  
end


# ==============================================================
# Solve for optimal drilling and production
# ==============================================================
function optpath(m::EOO_AKS, q0::Vector{Float64}, decreaseflag::Bool, 
    td0::Int64, x0::Vector{Float64}, mu0guess::Vector{Float64} = fill(2.0,NN(m)), 
    pvecg::Vector{Float64} = pguess(m))
    # INPUTS
    # m: an instance of EOO_AKS
    # q0: (Nx1) production in first period
    # decreaseflag: 0/1 flag. 0 = demand constant, 1 = demand decreases
    # td0: initial time period (positions the initial demand curve)
    # x0: (Nx1) initial capacity available to be drilled, mmbbl/d
    # mu0guess: (optional) (Nx1) guess of initial shadow value in $ per mmbbl/d
    # pvecg: (optional) initial price vector guess

    # OUTPUTS
    # D: (Nx1) cumulative drilling capacity added in mmbbl/d
    # Q: (Nx1) cumulative production in billion bbl
    # pvec: (Tx1) annual time series of prices in $/bbl
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period
    # qvec (TxN): rate of production in each period, mmbbl/d
    # mu0: (Nx1) initial (t=0) shadow value in $ per mmbbl/d
    # muvec: (TxN) time series of current shadow values in $ per mmbbl/d
    # Te (Nx1) periods when drilling stops for each region
    # Note: for pvec, dvec, qvec, and muvec, the first row corresponds to period td0

    # First address possibility that reserves are zero
    if x0 == zeros(NN(m))
        D = zeros(NN(m))
        dvec = zeros(maxT(m), NN(m))
        Te = zeros(NN(m))
        mu0 = zeros(NN(m))
        muvec = zeros(maxT(m), NN(m))
        pvec, qvec = pqpath(m, dvec, q0, decreaseflag, td0)
        Q = sum(qvec) * DY(m) / 1000 / T(m)
    else
        # See if D < reserves at mu0 = 0 for all fields
        D = cumdrill(m, zeros(NN(m)), q0, decreaseflag, td0)
        if all(D .<= x0)
            mu0 = zeros(NN(m))
        elseif all(x -> x == 0. || isinf(x), mu0guess)
            # mu0guess is all zero or Inf; don't search
            mu0 = mu0guess
        else
            # Need to search for mu0
            # Define utility function for the complementarity problem, lx is log(mu0)
            function utilc(m::EOO_AKS, lx::Vector{Float64}, q0::Vector{Float64},
                decreaseflag::Bool, td0::Int64, x0::Vector{Float64}, pvecg::Vector{Float64})
                cumD = cumdrill(m, exp.(lx), q0, decreaseflag, td0, pvecg) 
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
            function utilfm(m::EOO_AKS, lx::Vector{Float64}, q0::Vector{Float64}, decreaseflag::Bool, 
                td0::Int64, x0::Vector{Float64}, nzind::BitVector, infind::BitVector, pvecg::Vector{Float64})
                lxall = ones(NN(m))*(-20)  # initialize with small values
                lxall[infind] .= Inf
                lxall[nzind] = lx
                Dall = utilc(m, lxall, q0, decreaseflag, td0, x0, pvecg)
                return Dall[nzind]
            end
            # Solve for mu0
            result = nlsolve(lx -> utilfm(m, lx, q0, decreaseflag, td0, x0, nzind, infind, pvecg), 
                log.(mu0guess[nzind]), ftol=stol(m), show_trace=true, iterations=maxit(m))
            mu0 = zeros(NN(m))
            mu0[mu0guess.==Inf] .= Inf
            mu0[nzind] = exp.(result.zero)
            # Replace small values with zero
            mu0[mu0 .< mutol(m)] .= 0.
        end
        # Get time series of prices, drilling, and extraction
        pvec, dvec, qvec = solvepq(m, mu0, q0, decreaseflag, td0)
        D = vec(sum(dvec, dims=1))
        Q = vec(sum(qvec, dims=1)) * DY(m) / 1000 / T(m)
        Tvec = 0:maxT(m)-1
        muvec = mu0' .* (1 .+ r(m)') .^ Tvec
        # Stop times
        Te = zeros(NN(m))
        for n = 1:NN(m)
            if D[n] == 0
                Te[n] = 0
            else
                nz = findall(dvec[:, n] .> 0)
                Te[n] = maximum(nz)
            end
        end
    end
    return D, Q, pvec, dvec, qvec, mu0, muvec, Te
end


# ==============================================================
# Solve for myopic path in which firms don't believe demand is declining
# ==============================================================
function myopicpath(m::EOO_AKS, q0::Vector{Float64}, 
    mu0guess::Vector{Float64} = fill(2.0,NN(m)), pvecg::Vector{Float64} = pguess(m))
    # INPUTS
    # m: an instance of EOO_AKS
    # q0: (Nx1) production in first period
    # mu0guess: (optional) (Nx1) guess of initial shadow value in $ per mmbbl/d
    # pvecg: (optional) initial price vector guess

    # OUTPUTS
    # D: (Nx1) cumulative drilling capacity added in mmbbl/d
    # Q: (Nx1) cumulative production in billion bbl
    # pvec: (Tx1) annual time series of prices in $/bbl
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period
    # qvec (TxN): rate of production in each period, mmbbl/d
    # muvec: (TxN) time series of current shadow values in $ per mmbbl/d

    # Initialize loop
    t = 1
    endflag = false
    pvec = zeros(maxT(m)); dvec = zeros(maxT(m), NN(m)); qvec = zeros(maxT(m), NN(m))
    D = zeros(NN(m)); Q = zeros(NN(m))
    muvec = zeros(maxT(m), NN(m))
    xt = x0(m)    # initial reserves
    q0t = q0      # initial production rate
    lastp = min(shiftstop(m)+1, maxT(m))   # last period of loop
    
    # Loop over times
    while t <= lastp && endflag == false
        # Solve model assuming demand is fixed at its current level
        _, _, pt, dt, qt, mut, _, _ = optpath(m, q0t, false, t, xt, mu0guess, pvecg)
        # Store results
        pvec[t] = pt[1]; dvec[t, :] = dt[1, :]; qvec[t, :] = qt[1, :]; 
        muvec[t, :] = mut'
        # Decrease reserves
        xt -= vec(dt[1, :])
        println("Current t and value of xt: ", [t xt'])
        mu0guess = mut      # update guess
        mu0guess[xt .<= 0.] .= Inf  # set mu0guess to Inf if reserves are exhausted
        # Improve mu0 guess by forecasting ahead to the next iteration
        guessratio = ones(NN(m))
        if t >= 2
            for i = 1:NN(m)
                if muvec[t-1, i] > 0. && muvec[t, i] < Inf
                    guessratio[i] = mut[i] / muvec[t-1, i]                  
                end
            end
        end
        mu0guess = mu0guess .* guessratio
        println("Current mut and mu0guess: ", [mut' mu0guess'])
        # update initial production rate for next period
        q0t = vec(qt[2, :])
        # update price vector guess for next period
        pvecg = pt
        # Increment time
        t += 1
        # Trigger endflag if reserves are exhausted
        if maximum(xt) <= 0.
            endflag = true
        end
    end

    # Final set of periods when demand is at steady state again (possibly zero)
    if lastp==maxT(m) || drem(m)==0. || maximum(xt) <= 0.  
        # nothing else to do; all remaining periods have zero q and p
    else
        _, _, pte, dte, qte, mute, _, _ = optpath(m, q0t, false, t, xt, mu0guess, pvecg)
        pvec[t:end] = pte[1:maxT(m)-t+1]
        dvec[t:end, :] = dte[1:maxT(m)-t+1, :]
        qvec[t:end, :] = qte[1:maxT(m)-t+1, :]
        muvec[t:end, :] = repeat(mute', outer=[length(t:maxT(m)), 1])
    end
    
    # Cumulative drilling and production
    D = vec(sum(dvec, dims=1))
    Q = vec(sum(qvec, dims=1)) * DY(m) / 1000 / T(m)
    return D, Q, pvec, dvec, qvec, muvec
end


# ==============================================================
# Compute PV and profits from solved model
# ==============================================================
function cumprofits(m::EOO_AKS, pvec::Vector{Float64}, 
    dvec::Matrix{Float64}, qvec::Matrix{Float64})
    # INPUTS
    # m: an instance of EOO_AKS
    # pvec: (Tx1) annual time series of prices in $/bbl
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period
    # qvec (TxN): rate of production in each period, mmbbl/d

    # OUTPUTS
    # PVpi: (Nx1) PDV of profits for each region ($billion)
    # pivec: (TxN) time series of instantaneous profit flows ($billion / period)

    # Compute per-period profits (TxN)
    pivec = (pvec .* qvec * DY(m) / 1000 / T(m) 
        - α(m)' .* dvec - γ(m)' .* dvec.^2 / 2)

    # Compute PV
    Tvec = 0:maxT(m)-1
    Deltavec = δ(m)'.^Tvec
    Discpivec = pivec .* Deltavec
    PVpi = vec(sum(Discpivec, dims=1))
    
    return PVpi, pivec
end


# ==============================================================
# Compute first period drilling elasticity (no declining demand)
# ==============================================================
function drillelast(m::EOO_AKS, q0::Vector{Float64}, mu0guess=fill(2.0,NN(m)), 
    pvecg::Vector{Float64} = pguess(m))
    # INPUTS
    # m: an instance of EOO_AKS
    # q0: (Nx1) production in first period
    # mu0guess: (optional) (Nx1) guess of initial shadow value in $ per mmbbl/d
    # pvecg: (optional) initial price vector guess

    # OUTPUTS
    # elasttot: total drlg elasticity
    # elast: drlg elasticity for each region

    # Outcomes at actual first period demand
    _, _, pvec, dvec, _, _, _, _ = optpath(m, q0, false, 1, x0(m), mu0guess, pvecg)

    # Outcomes at shifted first period demand
    _, _, pvecs, dvecs, _, _, _, _ = optpath(m, q0, false, 2, x0(m), mu0guess, pvec)

    elast = ((dvecs[1, :] .- dvec[1, :]) ./ (pvecs[1] - pvec[1])  
        ./ (dvecs[1, :] .+ dvec[1, :]) .* (pvecs[1] + pvec[1]))
    elast = elast'

    elasttot = ((sum(dvecs[1, :]) - sum(dvec[1, :])) / (pvecs[1] - pvec[1]) ./ 
        (sum(dvecs[1, :]) + sum(dvec[1, :])) * (pvecs[1] + pvec[1]))

    return elasttot, elast
end


# ==============================================================
# Summarize quantities and profits
# ==============================================================
function sum_Q_pi(m::EOO_AKS, qvec::Matrix{Float64}, pvec::Vector{Float64}, 
    dvec::Matrix{Float64})
    # INPUTS
    # m: an instance of EOO
    # qvec: (TxN) time series of production rates in mmbbl/d
    # pvec: (Tx1) annual time series of prices in $/bbl
    # dvec (TxN): rate of drilling (capacity investment) in each period, mmbbl/d per period

    # OUTPUTS
    # QT: total production in billion bbl
    # PVQT: total discounted production in billion bbl
    # PVQ: discounted production by region, billion bbl
    # PVpi: PV of profits by region, $billion
    # elasttot: total supply elasticity
    # elast: supply elasticity for each region
    
    # Summaries of production
    QT = sum(qvec) * DY(m) / 1000 / T(m)     # total production in billion bbl
    # Discounted production by region, billion bbl
    PVQ = qvec' * Ediscvec(m) * DY(m) / 1000 / T(m)
    PVQT = sum(PVQ)                          # total discounted production in billion bbl

    # PV of profits and per period profits, by region
    PVpi, pivec = cumprofits(m, pvec, dvec, qvec)
    
    return QT, PVQT, PVQ, PVpi, pivec
end


# ==============================================================
# Return diff between grown q0 and second period production
# ==============================================================
function q0diff(m::EOO_AKS, q0in::Vector{Float64}, γguess::Vector{Float64},
    mu0guess::Vector{Float64}=fill(2.0,NN(m)))
    # INPUTS
    # m: an instance of EOO_N
    # q0in: (Nx1) actual first period production in mmbbl/d
    # γguess: (Nx1) guess of γ parameter vector, in $/bbl per mmbbl/d
    # mu0guess (optional): (Nx1) guess of initial shadow value in $ per mmbbl/d

    # OUTPUTS
    # diff: q0in minus first period production

    # Instantiate model with guessed γ
    mγ_temp = EOO_N(αin=α(m.eoo_n), γin=γguess, x0in=x0(m.eoo_n), delast=delast(m), 
        Pref=Pref(m), dref=dref(m), timetozero=timetozero(m), shiftdelin=shiftdel(m.eoo_n), 
        drem=drem(m), dgr_ann=dgr_ann(m), r_ann=r_ann(m), Er_ann=Er_ann(m), mp=mp(m), 
        maxY=maxY(m), DY=DY(m), T=T(m), stol=stol(m), mutol=mutol(m), maxit=maxit(m))
    mγ = EOO_AKS(eoo_n=mγ_temp, lambday=lambday(m), q0=q0in, pguess=pguess(m))

    # Simulate model and obtain production
    _, _, _, _, qvec, _, _, _ = optpath(mγ, q0in, false, 1, x0(mγ), mu0guess)

    # Grow q0in and take difference
    q0in2 = q0in .* exp.(dgr(m) .* min.(1, dgrt(m)))

    diff = q0in2 - vec(qvec[2, :])
    return diff
end


# ==============================================================
# Function to run AKS models given input parameter tuple
# ==============================================================
function runAKSmodel(pin::NamedTuple, unant::Bool = true)
    # INPUTS
    # pin: tuple of input parameters
    # unant: flag for whether to run unanticipated demand decline

    # OUTPUTS
    # Oi: tuple of baseline outputs
    # O: tuple of anticipated demand decline outputs
    # Om: tuple of myopic demand decline outputs
    # Op: tuple of percent changes outputs

    # Run non-AKS model to get guess of price vector and shadow value
    m_temp = EOO_N(αin=pin.alpha, γin=pin.gamma, x0in=pin.x0, 
        delast=pin.delast, Pref=pin.Pref, dref=pin.dref, 
        dgr_ann=pin.dgr_ann, timetozero=pin.timetozero, 
        r_ann=pin.r_ann, Er_ann=pin.Er_ann, mp=pin.mp, T=pin.T,
        shiftdelin=pin.shiftdelin, drem=pin.drem, maxY=pin.maxY, mutol=pin.mutol)
    _, pveci_temp, _, mu0i_temp, _, _ = optpath(m_temp, false, 0., x0(m_temp))
    # Now run AKS model
    m = EOO_AKS(eoo_n=m_temp, lambday=pin.lambday, q0=pin.q0, pguess=pveci_temp)
    Di, Qi, pveci, dveci, qveci, mu0i, muveci, Tei = 
        optpath(m, q0(m), false, 1, x0(m), mu0i_temp .* RRC(m))
    mci_bbl = mc(m, vec(dveci[1,:])) ./ RRC(m);             # initial MC in $/bbl
    mu0i_bbl = mu0i ./ RRC(m);                              # initial shadow val in $/bbl
    QTi, PVQTi, PVQi, PVpii, piveci = sum_Q_pi(m, qveci, pveci, dveci)
    # Demand decline
    Dd, Qd, pvecd, dvecd, qvecd, mu0d, muvecd, Ted = optpath(m, q0(m), true, 1, x0(m))
    QTd, PVQTd, PVQd, PVpid, pivecd = sum_Q_pi(m, qvecd, pvecd, dvecd)
    # Myopic
    if unant==true
        Dm, Qm, pvecm, dvecm, qvecm, muvecm = myopicpath(m, q0(m), mu0i, pveci);
        QTm, PVQTm, PVQm, PVpim, pivecm = sum_Q_pi(m, qvecm, pvecm, dvecm)
    else
        Dm = ones(NN(m)); Qm = ones(NN(m)); dvecm = ones(maxT(m),NN(m)); pvecm = ones(maxT(m),NN(m)); 
        qvecm = ones(maxT(m),NN(m)); muvecm = ones(maxT(m),NN(m));
        QTm = 1; PVQTm = 1; PVQm = ones(NN(m)); 
        PVpim = ones(NN(m)); pivecm = ones(maxT(m),NN(m));
    end
    # Percent changes
    PctIncDdDiTot = sum(Dd-Di) / sum(Di) * 100; 
    PctIncQdQiTot = (QTd-QTi) / QTi * 100; PctIncPVQdQiTot = (PVQTd-PVQTi) / PVQTi * 100
    if unant==true
        PctIncQmQiTot = (QTm-QTi) / QTi * 100; PctIncPVQmQiTot = (PVQTm-PVQTi) / PVQTi * 100
        PctIncDdDmTot = sum(Dd-Dm) / sum(Dm) * 100; 
        PctIncQdQmTot = (QTd-QTm) / QTm * 100; PctIncPVQdQmTot = (PVQTd-PVQTm) / PVQTm * 100
    else
        PctIncQmQiTot = 1; PctIncPVQmQiTot = 1; PctIncDdDmTot = 1; 
        PctIncQdQmTot = 1; PctIncPVQdQmTot = 1
    end
    # Create output touples
    Oi = (D = Di, Q = Qi, pvec = pveci, dvec = dveci, qvec = qveci, 
        mu0 = mu0i, muvec = muveci, Te = Tei,
        mc_bbl = mci_bbl, mu0_bbl = mu0i_bbl, QT = QTi, 
        PVQT = PVQTi, PVQ = PVQi, PVpi = PVpii, pivec = piveci)
    Od = (D = Dd, Q = Qd, pvec = pvecd, dvec = dvecd, qvec = qvecd,
        mu0 = mu0d, muvec = muvecd, Te = Ted, QT = QTd, 
        PVQT = PVQTd, PVQ = PVQd, PVpi = PVpid, pivec = pivecd)
    Om = (D = Dm, Q = Qm, pvec = pvecm, dvec = dvecm, qvec = qvecm,
        muvec = muvecm, QT = QTm, PVQT = PVQTm, PVQ = PVQm, 
        PVpi = PVpim, pivec = pivecm)
    Op = (PctIncDdDiTot = PctIncDdDiTot, PctIncDdDmTot = PctIncDdDmTot, 
        PctIncQdQmTot = PctIncQdQmTot, PctIncPVQdQmTot = PctIncPVQdQmTot, 
        PctIncQdQiTot = PctIncQdQiTot, PctIncPVQdQiTot = PctIncPVQdQiTot,
        PctIncQmQiTot = PctIncQmQiTot, PctIncPVQmQiTot = PctIncPVQmQiTot)
    return Oi, Od, Om, Op
end