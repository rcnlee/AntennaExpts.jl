using AntennaExpts, YagiAntProblem
using ExprSearch, DerivTreeVis
using ExprSearch.GP
using RLESUtils
import RLESTypes.SymbolTable

function test()
    problem = YagiAnt()
    grammar = get_grammar(problem)
    mda = min_depth_actions(grammar) 
    ind = rand(grammar, mda, 10)
    derivtreevis(ind.derivtree, "treevis")
    expr = get_expr(ind.derivtree)
    fitness = get_fitness(problem, ind.derivtree, SymbolTable()) 
end
