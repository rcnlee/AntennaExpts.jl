@everywhere using AntennaExpts
@everywhere using PIFAProblem
using PIFA_GE
using ExprSearch
using ExprSearch.GE

function ExprSearch.GE.evaluate!(p::GEESParams, pop::GEPopulation, result::GEESResult, 
    problem::PIFAAnt, default_expr)
    #first convert all derivtrees to exprs
    for ind in pop
        if isnull(ind.fitness)
            try
                ind.expr = get_expr(ind.derivtree)
            catch e
                if !isa(e, IncompleteException)
                    rethrow(e)
                end
                ind.expr = default_expr
                ind.fitness = Nullable{Float64}(realmax(Float64))
            end
        end
    end
    asyncdriver(pop, result)
end

function asyncdriver(pop::GEPopulation, result::GEESResult)
    np = nprocs()  # determine the number of processes available
    n = length(pop)
    i = 1
    function nextidx() 
        while i <= n 
            idx = i
            i += 1
            if isnull(pop[idx].fitness)
                return idx
            end
        end
        n+1
    end
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    #println("myid=", myid())
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        #println("idx=", idx, ", p=", p)
                        fitness = remotecall_fetch(worker_expr_fitness, p, pop[idx].expr)
                        pop[idx].fitness = Nullable{Float64}(fitness)
                        result.totalevals += 1
                        if fitness < result.fitness
                            result.fitness = fitness
                            copy!(result.tree, pop[idx].derivtree)
                            result.best_at_eval = result.totalevals
                            result.expr = pop[idx].expr
                        end
                    end
                end
            end
        end
    end
end

@everywhere problem = YagiAnt()
@everywhere worker_expr_fitness(expr) = begin 
    #println("myid=", myid()) 
    expr_fitness(problem, expr)
end

pifa_ge(; seed=1, pop_size=20000, tournament_size=4000, iterations=20, top_keep=0.0)
