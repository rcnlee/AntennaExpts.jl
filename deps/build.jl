pkgs = Pkg.installed()
!haskey(pkgs, "RLESUtils") && Pkg.clone("https://github.com/rcnlee/RLESUtils.jl.git", "RLESUtils")
!haskey(pkgs, "NECPP") && Pkg.clone("https://github.com/rcnlee/NECPP.jl.git", "NECPP")
