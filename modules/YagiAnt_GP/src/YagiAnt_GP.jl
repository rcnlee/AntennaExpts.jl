"""
Yagi-like antenna fitting inside a Bott's dot
"""
module YagiAnt_GP

export yagiant_gp

using ExprSearch, YagiAntProblem, DerivTreeVis
using ExprSearch.GP
using RLESUtils, LogSystems, Loggers

const DIR = dirname(@__FILE__)
const RESULTDIR = joinpath(DIR, "..", "..", "..", "results")

"""
Example usage:
yagiant_gp()
"""
function yagiant_gp(; outdir::AbstractString=joinpath(RESULTDIR, "YagiAnt_GP"),
                     seed=1,
                     logfileroot::AbstractString="yagiant_gp_log",

                     pop_size::Int64=1000,
                     maxdepth::Int64=10,
                     iterations::Int64=10,
                     tournament_size::Int64=20,
                     top_keep::Float64=0.01,
                     crossover_frac::Float64=0.4,
                     mutate_frac::Float64=0.2,
                     rand_frac::Float64=0.2,
                     default_code::Any=0.0,

                     departure_angle::Float64=10.0,
                     vis::Bool=true)

    srand(seed)
    mkpath(outdir)

    logsys = GP.logsystem()
    send_to!(STDOUT, logsys, ["verbose1", "current_best_print", "result"])
    logs = TaggedDFLogger()
    send_to!(logs, logsys, ["code", "computeinfo", "current_best", "elapsed_cpu_s", "fitness",
        "fitness5", "parameters", "result"])

    problem = YagiAnt(departure_angle)
  
    gp_params = GPESParams(pop_size, maxdepth, iterations, tournament_size, top_keep,
        crossover_frac, mutate_frac, rand_frac, default_code, logsys)
  
    result = exprsearch(gp_params, problem)

    outfile = joinpath(outdir, "$(logfileroot).txt")
    save_log(outfile, logs)

    if vis
        #visualize the deriv tree
        derivtreevis(get_derivtree(result), joinpath(outdir, "$(logfileroot)_derivtreevis"))

        #output the nec file
        f = open("result.nec", "w")
        (gain, swr) = nec_run(problem, result.expr, f) do nec
            gain = nec_gain(nec, 0, 0, 0)
            imp_re = nec_impedance_real(nec, 0)
            imp_im = nec_impedance_imag(nec, 0)
            Z = Complex(imp_re, imp_im)
            Z0 = Complex(50.0, 0.0)
            swr = vswr(Z, Z0)
            (gain, swr) 
        end
        close(f)
        @show (gain, swr)
    end
    @show result.expr
    
    return result
end

end #module
