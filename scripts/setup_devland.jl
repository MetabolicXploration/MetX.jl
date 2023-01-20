import Pkg

## ------------------------------------------------------------------
const MET_X_PKGS = [
    "MetX", 
    "MetXBase", 
    "MetXOptim", 
    "MetXGrids", "MetXMC", "MetXEP", 
    "MetXPlots", 
    "MetXNetHub"
]
const BASE_URL = "https://github.com/MetabolicXploration"

## ------------------------------------------------------------------
# add registry
let
    println("\n\n")
    println("="^60)
    println("DOWNLOADING")
    
    url = joinpath(BASE_URL, "MetX_Registry_jl")
    Pkg.Registry.add(Pkg.RegistrySpec(;url))
end

## ------------------------------------------------------------------
# try download to Pkg.pkgdir
let
    println("\n\n")
    println("="^60)
    println("DOWNLOADING")
    
    for pkgname in MET_X_PKGS
        
        pkgdir = joinpath(Pkg.devdir(), pkgname)
        url = joinpath(BASE_URL, string(pkgname, ".jl"))
        
        println()
        println("."^60)
        println()
        @info("Doing", pkgname, url, pkgdir)
        println()
        
        if isdir(pkgdir) 
            @info("Package already downloaded")
            continue
        end

        # Download
        run(`git clone $(url) $(pkgdir)`)
    end
end

## ------------------------------------------------------------------
# try to pull
let
    
    println("\n\n")
    println("="^60)
    println("PULLING")

    for pkgname in MET_X_PKGS
        
        pkgdir = joinpath(Pkg.devdir(), pkgname)
        url = joinpath(BASE_URL, string(pkgname, ".jl"))
        
        println()
        println("."^60)
        println()
        @info("Doing", pkgname, url, pkgdir)
        println()
        
        mkpath(pkgdir)
        cd(pkgdir) do

            if !success(`git diff --quiet`)
                @warn("Package repo is dirty")
                println()
                run(`git -C $(pkgdir) status`)
                
                println()
                @info("Trying pull --ff-only")
                println()

                if success(`git -C $(pkgdir) pull --ff-only`)
                    println()
                    @info("pull --ff-only OK")
                else
                    println()
                    @warn("pull --ff-only Failed")
                end
            else
                # Pull
                if success(`git -C $(pkgdir) pull`)
                    println()
                    @info("pull OK")
                else
                    println()
                    @info("pull Failed")
                end
            end
        end

    end
end

## ------------------------------------------------------------------
# Dev in home pkg
let
    println("\n\n")
    println("="^60)
    println("DEV IN HOME")

    Pkg.activate()
    for pkgname in MET_X_PKGS

        println()
        println("."^60)
        println()
        @info("Doing", pkgname)
        println()
        
        Pkg.develop(pkgname)
        
        println()
    end
end

## ------------------------------------------------------------------
# Dev cross-dependencies
function _get_deps(projtoml) 
    tomldict = Pkg.TOML.parsefile(projtoml)
    return get(tomldict, "deps", Dict()) |> keys |> collect
end

let
    println("\n\n")
    println("="^60)
    println("DEV CROSS DEPENDENCIES")

    Pkg.activate()
    for pkgname in MET_X_PKGS

        pkgdir = joinpath(Pkg.devdir(), pkgname)
        projtoml = joinpath(pkgdir, "Project.toml")
        
        if !isfile(projtoml) 
            @warn("Project.Toml no found. Skipping", projtoml)
            continue
        end
        
        println()
        println("."^60)
        println()
        @info("Doing", pkgname, projtoml)
        println()
        
        Pkg.activate(projtoml)
        println()
        Pkg.status()
        println()

        deps = _get_deps(projtoml)

        for pkgname1 in MET_X_PKGS
            (pkgname1 == pkgname) && continue
            (pkgname1 in deps) || continue

            Pkg.develop(pkgname1)
            println()
        end

        # instantiate
        println()
        Pkg.instantiate()
        println()
    end

end

## ------------------------------------------------------------------
println("\n\n")
println("="^60)
println("DONE!!! ENJOY")
println("\n\n")


