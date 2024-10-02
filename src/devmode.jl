const METX_PKGS = [
    "MetXBase", "MetXGEMs", "MetXOptim",
    "MetXEP", "MetXMC", "MetXNetHub", 
    "MetXCultureHub"
]

# Make all MetX packages be in dev mode
function devmode()
    for pkgname in METX_PKGS
        Pkg.develop(pkgname)
    end
end