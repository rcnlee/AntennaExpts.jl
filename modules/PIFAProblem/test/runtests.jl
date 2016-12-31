using AntennaExpts, PIFAProblem
using ExprSearch, DerivTreeVis
using ExprSearch.GE
using RLESUtils
import RLESTypes.SymbolTable
using NECPP

srand(1)
function test(N)
    problem = PIFAAnt()
    grammar = get_grammar(problem)
    mda = min_depth_actions(grammar) 
    for i = 1:N
        ind = rand(grammar, mda, 10)
        #derivtreevis(ind.derivtree, "treevis")
        expr = get_expr(ind.derivtree)
        f = open("test.nec","w")
        (gain, Z) = nec_run(problem, expr, f) do nec
            gain = nec_gain(nec, 0, 0, 0)
            imp_re = nec_impedance_real(nec, 0)
            imp_im = nec_impedance_imag(nec, 0)
            Z = Complex(imp_re, imp_im)
            (gain, Z) 
        end
        close(f)
        Z0 = Complex(50.0, 0.0)
        v = vswr(Z, Z0)
        @show gain
        @show Z
        @show v
        v < 0.0 && break
    end
end

