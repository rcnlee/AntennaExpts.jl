module AntennaExpts

const MODULEDIR = joinpath(dirname(@__FILE__), "..", "modules")

using RLESUtils.ModLoader

load_to_path(MODULEDIR)
const PKGS = readdir(MODULEDIR)

"""
Test an individual submodule
"""
function test(pkgs::AbstractString...; coverage::Bool=false)
  cd(() -> Pkg.Entry.test(AbstractString[pkgs...]; coverage=coverage), MODULEDIR)
end

"""
Test all submodules in modules folder.  Don't stop on error.
"""
function testall()
    for pkg in PKGS
        try
            test(pkg)
        catch
            println("Error in $pkg")
        end
    end
end

end # module
