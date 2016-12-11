using AntennaExpts, YagiAntProblem, YagiAnt_GP
using ExprSearch
using ExprSearch.GP

function ExprSearch.GP.evaluate!(p::GPESParams, pop::GPPopulation, result::GPESResult, 
    problem::YagiAnt, default_expr)
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

function asyncdriver(pop::GPPopulation, result::GPESResult)
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

@everywhere using YagiAntProblem
@everywhere problem = YagiAnt()
@everywhere worker_expr_fitness(expr) = begin 
    #println("myid=", myid()) 
    expr_fitness(problem, expr)
end

yagiant_gp(; seed=2, pop_size=10000, tournament_size=1000, iterations=20)
