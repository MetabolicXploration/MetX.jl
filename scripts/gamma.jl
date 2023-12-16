## ----------------------------------------------------
@time begin 
    using MetXBase
    using MetXEP
    using MetXEP: gamma_gd!
    using MetXNetHub
    using MetXOptim
    using MetXEP
    using Gurobi
    using Optim
end

## ------------------------------------------------------------------
LIN_SOLVER = Gurobi.Optimizer

# ------------------------------------------------------------------
# gamma grad desc
let
    ## ---------------------------------------------
    model_id = "ecoli_core"
    net0 = pull_net(model_id)
    global lep = box(net0, LIN_SOLVER; verbose = false)
    global biom_id = extras(lep, "BIOM") 

    ## ---------------------------------------------
    # EP MODEL
    global epm = FluxEPModelT0(lep)
    config!(epm; epsconv = 1e-4, verbose = true)
    converge!(epm)

    # ## ---------------------------------------------
    target_ids = [biom_id]
    # target_values = [0.002]
    target_values = [0.0006]
    obj_fun = (x...) -> var(epm, target_ids)
    gdmodel = gamma_gd!(obj_fun, epm, target_ids, target_values;
        maxΔx = [1e6],
        minΔx = [5e5],
        gdth = 1e-3,
        maxiter = 1500,
        verbose = true
    )
    @show var(epm, target_ids)

    ## ---------------------------------------------
    @assert isapprox(obj_fun(), target_values; rtol = 1e-2)
end