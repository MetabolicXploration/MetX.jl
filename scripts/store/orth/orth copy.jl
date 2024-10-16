@time begin
    using LinearAlgebra
    using Plots
    using MetXEP
    using MetXNetHub
    using MetXOptim
    using MetXBase
    using MetXBase: mgrscho
    using MetXOptim: GLPK
    using Plots
end

## ---------------------------------------------
let
    model_id = "ecoli_core"
    global net0 = pull_net(model_id)
    biom_id = extras(net0, "BIOM")
    glc_id = extras(net0, "EX_GLC")
    
    global net = fva_strip(net0, GLPK.Optimizer)
    global epm = FluxEPModelT0(net)
    config!(epm, :verbose, true)
    
    betas = 10.0.^range(-5, 7.0; length = 50)
    Ss = Float64[]
    bioms = Float64[]
    for b in betas

        @info("Hello", b)
        
        beta!(epm, biom_id, b)
        converge!(epm)
        push!(bioms, mean(epm, biom_id))
        
        # S = entropy(epm, net)
        # push!(Ss, S)
    
    end

    # plot(betas, Ss; label = "")
    plot(log10.(betas), bioms; label = "")
end

## ---------------------------------------------
# function entropy(epm, net)
#     idxf, idxd, C, y, basis = _echelonize2(net.S, net.b)
#     # basis * vi = v -> basis
#     ort_basis = mgrscho(basis)
#     Cov = basis * epm.Σi * basis'
#     normΣ = ort_basis' * Cov * ort_basis
#     normΣ .= 0.5 * (normΣ + normΣ') # Correct Symetry
    
#     N = size(normΣ,1);
#     # this is just the entropy of a multivariate gaussian
#     L = cholesky(normΣ).L;
#     S = sum(log.(diag(L))) + 0.5*N*log.(2*pi*exp(1));
#     return S 
# end


