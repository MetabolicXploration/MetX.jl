. # --------------------------------------------------------------------
# globals
MET_X_PKG_NAMES = [
    "MetXBase", 
    "MetXGEMs", 
    "MetXOptim", 
    "MetXGrids", "MetXMC", "MetXEP", 
    "MetXPlots", "RunTestsEnv",
    "MetXNetHub", "MetXCultureHub", 
    "MetX", 
]
METX_ORG_URL = "https://github.com/MetabolicXploration"

# --------------------------------------------------------------------
# Utils
function _foreach_dep(f::Function)
    deps = Pkg.dependencies()
    for dep in values(deps)
        f(dep) == true && return 
    end
end

# --------------------------------------------------------------------
nothing