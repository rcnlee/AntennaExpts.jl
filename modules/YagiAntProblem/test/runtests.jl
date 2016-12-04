using AntennaExpts, YagiAntProblem
using ExprSearch, DerivTreeVis
using ExprSearch.GP
using RLESUtils
import RLESTypes.SymbolTable
using NECPP

function test()
    problem = YagiAnt()
    grammar = get_grammar(problem)
    mda = min_depth_actions(grammar) 
    ind = rand(grammar, mda, 10)
    derivtreevis(ind.derivtree, "treevis")
    expr = get_expr(ind.derivtree)
    f = open("test.nec","w")
    fitness = nec_run(problem, expr, f) do nec
        #for now... look at impedance
        result_index = 0
        z = Complex(nec_impedance_real(nec,result_index), nec_impedance_imag(nec,result_index))
        abs(z) #magnitude 
    end
    close(f)
    fitness
end

