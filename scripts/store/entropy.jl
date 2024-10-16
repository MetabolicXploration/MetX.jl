@time begin
    using LinearAlgebra
    using LinearAlgebra: eigen, Diagonal, I, det, rank
    using Plots
    using MetXEP
    using MetXEP: entropy
    using MetXNetHub
    using MetXOptim
    using MetXOptim: GLPK
    using MetXBase
    using MetXBase: mgrscho, nearPD, nearPD!
    using Plots
end

## ---------------------------------------------
let

    for _ in 1:50
        A = rand(100, 100)
        A .= 0.5 * (A + A')
        @assert !isposdef(A)
        nearPD!(A)
        @assert isposdef(A)
    end

end
## ---------------------------------------------

function _round!(epm; digits = 10)
    for f in fieldnames(typeof(epm))
        v = getfield(epm, f)
        typeof(v) <: Array || continue
        eltype(v) <: Number || continue
        @inbounds for i in eachindex(v)
            v[i] = round(v[i]; digits)
        end
    end
    return epm
end


## ---------------------------------------------
# MetNet
@time let
    global model_id = "iJO1366"
    global net0 = pull_net(model_id)
    global biom_id = extras(net0, "BIOM")
    global glc_id = extras(net0, "EX_GLC")
    global net0 = fva_strip(net0, GLPK.Optimizer; 
        protect_obj = true, verbose = true
    )
    todel = ["R_RBFSb" ,"R_PMDPHT" ,"R_GTPCII2" ,"R_DHPPDA2" ,"R_DB4PS" ,"R_RBFSa" ,"R_DPCOAK" ,"R_PPNCL2" ,"R_ASP1DC" ,"R_DPR" ,"R_PTPATi" ,"R_PANTS" ,"R_MOHMT" ,"R_PPCDC" ,"R_PNTK"]
    # empty_rxn!(net0, todel)
    # global net0 = emptyless_model(net0)

    global opm = fba(net0, GLPK.Optimizer)
    @show solution(opm, biom_id)
end

## ---------------------------------------------
# EP
let
    global epm0 = FluxEPModelT0(net0)
    config!(epm0; 
        verbose = true,
        epsconv = 1e-6, 
        maxvar = 1e15,
        minvar = 1e-15,
        damp = 0.99,
    )
    return nothing
end

## ---------------------------------------------
let
    r = abs.(net0.lb .- net0.ub)
    idxs = sortperm(r)
    println.(net0.rxns[idxs], " - ", r[idxs], " | ", net0.lb[idxs])
    nothing
end

## ---------------------------------------------
let
    global epm = FluxEPModelT0(net0)

    global errs = Float64[]
    global biomsv = Float64[]
    oniter = (epm) -> begin
        push!(errs, state(epm, :max_err))
        push!(biomsv, mean(epm, biom_id))
    end

    config!(epm; 
        verbose = true,
        epsconv = 1e-5, 
        maxiter = 10000, 
        # maxvar = 1e10,
        # minvar = 1e-10,
        damp = 0.94,
    )
    beta!(epm, biom_id, 50.0)
    converge!(epm; oniter)

    plot(biomsv)

end

## ---------------------------------------------
plot()
## ---------------------------------------------
## ---------------------------------------------
## ---------------------------------------------
## ---------------------------------------------
## ---------------------------------------------
## ---------------------------------------------
## ---------------------------------------------
@time let

    global net = deepcopy(net0)
    global epm = deepcopy(epm0)
    Biomass_max = ub(net, biom_id)    
    
    global betas = 10.0.^range(-5, 6.0; length = 50)
    global Ss = Float64[]
    global bioms = Float64[]
    global flxs = zeros(length(betas), size(net, 2)) # sample, flx
    oniter = (epm) -> begin
        # nearPD!(epm.Σi, 1e-10)
        # nearPD!(epm.Σd, 1e-10)
        return nothing
    end
    for (i, b) in enumerate(betas)
        
        beta!(epm, biom_id, b)
        converge!(epm; oniter)

        biom = mean(epm, biom_id)
        push!(bioms, biom)
        flxs[i, :] .= mean(epm)
        
        S = 1
        S = entropy(epm)
        push!(Ss, S)

        @info("Done", 
            beta = beta(epm, biom_id),
            status = state(epm, :status), 
            iter = state(epm, :iter), 
            Entropy = S,
            Biomass = biom,
            Biomass_max = ub(net, biom_id)
        )
    
    end
    return nothing
end

## ---------------------------------------------
function _normalize(v)
    v = v .- minimum(v)
    v ./= maximum(abs, v)
    return v
end

let
    p = plot(; 
        title = string(model_id, " ", size(net0)),
        legend = :left, 
        xlabel = "log(beta)",
        ylabel = "\$Normalized~ \\frac{v - min(v)}{max(abs(v))}\$",
    )
    
    plot!(p, log10.(betas), _normalize.(eachcol(flxs)); label = "", 
        c = :gray, alpha = 0.5, lw = 2
    )
    plot!(p, log10.(betas), _normalize(flxs[:, 1]),
        c = :gray, alpha = 0.5, lw = 2,
        label = "All Reactions"
    )
    plot!(p, log10.(betas), _normalize(Ss); 
        label = "Entropy", lw = 3, 
        c = :black
    )
    plot!(p, log10.(betas), _normalize(bioms); 
        label = "Biomass", lw = 3,
        c = :blue
    )
    p
end