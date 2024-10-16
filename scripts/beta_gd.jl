## ----------------------------------------------------
@time begin 
    using MetXBase
    using MetXEP
    using MetXEP: average_gd!
    using MetXNetHub
    using MetXOptim
    using MetXEP
    using Gurobi
    using Optim
end

## ---------------------------------------------
let

    ## ---------------------------------------------
    model_id = "ecoli_core"
    global net0 = pull_net(model_id)
    global lep = fva_strip(net0, Gurobi.Optimizer; nths = 2, verbose = true)
    biom_id = extras(lep, "BIOM") # R_BIOMASS_Ecoli
    glc_id = extras(lep, "EX_GLC") # R_EX_glc__D_e
    lac_id = "EX_lac__D_e"

    ## ---------------------------------------------
    # target_ids = [biom_id, glc_id, lac_id]
    target_ids = [biom_id]
    # target_value = [0.2, -9.0, 0.1]
    target_value = [0.2]
    beta0 = zeros(length(target_ids)) .+ 1e2

    ## ---------------------------------------------
    # EP MODEL
    global epm = FluxEPModelT0(lep)
    config!(epm; epsconv = 1e-4, verbose = true)
    converge!(epm)
    config!(epm; verbose = false)

    ## ---------------------------------------------
    # GD
    solver = GradientDescent()
    # solver = BFGS()
    # solver = LBFGS()
    # solver = ConjugateGradient()
    solver = Newton()
    # solver = NelderMead()
    opt = Optim.Options(show_trace = false, iterations=10000)
    result = optimize(beta0, solver, opt) do betai
        # Up betas
        beta!(epm, target_ids, betai)
        # converge
        converge!(epm)
        # lost function
        aves = mean(epm, target_ids)
        err = mean((aves .- target_value).^2)
        @show (betai, aves, err)
        return err
    end

    @show result.minimizer
    @show result.minimum

    return nothing
end

## ---------------------------------------------
let 
    ## ---------------------------------------------
    model_id = "ecoli_core"
    global net0 = pull_net(model_id)
    global lep = fva_strip(net0, Gurobi.Optimizer; nths = 2, verbose = true)
    biom_id = extras(lep, "BIOM") # R_BIOMASS_Ecoli
    glc_id = extras(lep, "EX_GLC") # R_EX_glc__D_e
    lac_id = "EX_lac__D_e"
    
    ## ---------------------------------------------
    # EP MODEL
    global epm = FluxEPModelT0(lep)
    config!(epm; epsconv = 1e-4, verbose = true)
    converge!(epm)
    
    ## ---------------------------------------------
    target_ids = [biom_id, glc_id, lac_id]
    target_value = [0.2, -9.0, 0.1]
    maxΔx = [1e3, 1e3, 1e4]
    gdmodel = average_gd!(epm, target_ids, target_value;
        maxΔx,
        gdth = 1e-2,
        maxiter = 500,
        verbose = true
    )
    
    ## ---------------------------------------------
    # @assert isapprox(up_fun(gdmodel), target; atol = 1e-3)
    @test isapprox(mean(epm, target_ids), target_value; rtol = 1e-2)
    println()
end

# Amoebacrew

## ---------------------------------------------
let
    TH_TESTS_LINSOLVER = Clp.Optimizer
    
    verbose = true

    model_id = "iJR904"
    net0 = pull_net(model_id)
    lep0 = lepmodel(net0)
    @time lep2 = fva_strip(lep0, TH_TESTS_LINSOLVER; nths = 1, verbose)
    @time lep2 = fva_strip(lep0, TH_TESTS_LINSOLVER; nths = 2, verbose)
    @time lep1 = fva_strip(lep0, TH_TESTS_LINSOLVER; nths = 3, verbose)
    @time lep1 = fva_strip(lep0, TH_TESTS_LINSOLVER; nths = 4, verbose)
end