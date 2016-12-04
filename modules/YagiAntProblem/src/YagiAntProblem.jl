"""
Yagi-like antenna fitting inside a Bott's dot
"""
module YagiAntProblem

export YagiAnt, create_grammar, get_grammar, get_fitness, 
    OutOfBoundsException, nec_run

using ExprSearch, DerivationTrees
import ExprSearch: ExprProblem, get_fitness, get_grammar
using RLESUtils, Interpreter, LogSystems
using NECPP

#all measurements in meters
const MAXLENGTH = 0.12
const MAXWIDTH = 0.11
const N_SEGMENTS = 11
const WIRE_RADIUS = 0.001 #meters
const FREQ_MHZ = 915.0 #MHz
const LAMBDA = 299.79 / FREQ_MHZ
const EXTYPE = 0 #0=voltage source

const HEIGHT = 0.02
const WXTOL = 0.02 #buffer space around wx
const TOL = 0.005 #min dist between two wires
const WIRE_LEN_TOL = 0.2 #antenna length tolerance from lambda/2
const WIRE_LEN_MIN = (1.0 - WIRE_LEN_TOL) * LAMBDA / 2
const WIRE_LEN_MAX = (1.0 + WIRE_LEN_TOL) * LAMBDA / 2 

const DIVCONST = 100.0

immutable OutOfBoundsException <: Exception end

type YagiAnt <: ExprProblem
    departure_angle::Float64 #degrees rel horizontal
    grammar::Grammar
    Z0::Complex
end

function YagiAnt(departure_angle::Float64=10.0, 
    Z0::Complex=Complex(50.0, 0.0))
    grammar = create_grammar()
    YagiAnt(departure_angle, grammar, Z0)
end

function create_grammar()
    @grammar grammar begin
        start = ant
        ant[make_block] = assign_exc + dir_el
        assign_exc = Expr(:(=), :lxwx, exc_el)
        exc_el = Expr(:call, :exc_wire, :nec, rn, rn)
        dir_el = dir_case2a | dir_case2b
        dir_case2a = Expr(:call, :dir_wire_2a, :nec, :lxwx, rn, rn, rn)
        dir_case2b = Expr(:call, :dir_wire_2b, :nec, :lxwx, rn, rn, rn)
        rn[divconst] = 0:100 #fraction of some length
    end
    grammar
end

###
#functions for grammar
make_block(lst::Array) = Expr(:block, lst...)

function exc_wire(nec::NecContext, rn_tl::Float64, rn_w::Float64)
    tlx = WIRE_LEN_MIN + rn_tl * (WIRE_LEN_MAX - WIRE_LEN_MIN)
    wxmin = WXTOL
    wxmax = MAXWIDTH
    wx = wxmin + rn_w * (wxmax - wxmin) 
    lx = (tlx - wx) / 2
    sym_bent_wire(nec, 1, 0.0, HEIGHT, WIRE_RADIUS, tlx, 
        wx, true)
    (lx, wx)
end

function dir_wire_2a(nec::NecContext, lxwx::Tuple, rn_tl::Float64,
    rn_w::Float64, rn_xpos::Float64)

    lx, wx = lxwx
    tld = WIRE_LEN_MIN + rn_tl * (WIRE_LEN_MAX - WIRE_LEN_MIN)

    wdmin = TOL 
    wdmax = wx - TOL 
    wd = wdmin + rn_w * (wdmax - wdmin) 
    ld = (tld - wd) / 2

    xposmin = ld + TOL
    xposmax = MAXLENGTH 
    xpos = xposmin + rn_xpos * (xposmax - xposmin)

    sym_bent_wire(nec, 4, xpos, HEIGHT, WIRE_RADIUS, tld, 
        wd, false)
    nothing
end
function dir_wire_2b(nec::NecContext, lxwx::Tuple, rn_tl::Float64,
    rn_w::Float64, rn_xpos::Float64)

    lx, wx = lxwx
    tld = WIRE_LEN_MIN + rn_tl * (WIRE_LEN_MAX - WIRE_LEN_MIN)

    wdmin = TOL 
    wdmax = wx - TOL 
    wd = wdmin + rn_w * (wdmax - wdmin) 
    ld = (tld - wd) / 2

    xposmin = TOL
    xposmax = MAXLENGTH - ld
    xpos = xposmin + rn_xpos * (xposmax - xposmin)

    sym_bent_wire(nec, 4, xpos, HEIGHT, WIRE_RADIUS, tld, 
        wd, true)
    nothing
end
function sym_bent_wire(nec::NecContext, tag_id_base::Int64, xpos::Float64, 
    height::Float64, wire_radius::Float64, total_len::Float64, 
    width::Float64, up::Bool)

    halfwidth = width / 2
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    zw1 = zw2 = height
    rad = wire_radius

    tag_id = tag_id_base 
    segment_count = N_SEGMENTS
    xw1 = xpos
    yw1 = +halfwidth 
    xw2 = xpos 
    yw2 = -halfwidth
    check_boundaries([xw1, xw2], 0.0, MAXLENGTH)
    check_boundaries([yw1, yw2], -MAXWIDTH/2, MAXWIDTH/2)
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad)) #base

    for y in [+halfwidth, -halfwidth]
        tag_id += 1 
        segment_count = N_SEGMENTS
        xw1 = xpos
        yw1 = y 
        xlen = (total_len - width)/2 
        xw2 = up ? xpos + xlen : xpos - xlen
        yw2 = y
        check_boundaries([xw1, xw2], 0.0, MAXLENGTH)
        check_boundaries([yw1, yw2], -MAXWIDTH/2, MAXWIDTH/2)
        handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
            xw2, yw2, zw2, rad, rdel, rrad)) 
    end
end

function check_boundaries(xs::Vector{Float64}, xmin::Float64, xmax::Float64)
    for x in xs
        if x < xmin || x > xmax
            throw(OutOfBoundsException())
        end
    end
end

divconst(x) = x / DIVCONST

const SYMTABLE = SymbolTable(
    :make_block => make_block,
    :exc_wire => exc_wire,
    :dir_wire_2a => dir_wire_2a,
    :dir_wire_2b => dir_wire_2b,
    :divconst => divconst
    )

function eval_expr(problem::YagiAnt, expr, nec) 
    SYMTABLE[:nec] = nec
    interpret(SYMTABLE, expr)
end

ExprSearch.get_grammar(problem::YagiAnt) = problem.grammar

function ExprSearch.get_fitness(problem::YagiAnt, derivtree::DerivationTree, 
    userargs::SymbolTable)

    expr = get_expr(derivtree)
    fitness = nec_run(expr) do nec
        gain = nec_gain(nec, 0, 0, 0)
        imp_re = nec_impedance_real(nec, 0)
        imp_im = nec_impedance_imag(nec, 0)
        Z = Complex(imp_re, imp_im)
        gain + vswr_metric(vswr(Z, problem.Z0))
    end
    fitness
end

#vswr_metric(x::Float64) = 10.0 / (1+e^-(2x-5))
function vswr_metric(x::Float64) 
    #passes through (0,0), (1.5,0.5), (3,10)
    if 0.0 <= x <= 1.5
        return x / 3
    elseif 1.5 < x <= 3.0
        return 6.3333x - 9.0
    elseif 3.0 < x
        return 10.0
    else
        throw(DomainError("vswr cannot be negative!"))
    end
end

nec_run(f::Function, problem::YagiAnt, expr, io::IO) = nec_run(f, problem, expr, Nullable{IO}(io))
function nec_run(f::Function, problem::YagiAnt, expr, io::Nullable{IO}=Nullable{IO}())
    if !isnull(io)
        logsys = NECPP.logsystem()
        send_to!(get(io), logsys, "nec_cards", first)
        nec = nec_create(logsys)
    else
        nec = nec_create()
    end
    eval_expr(problem, expr, nec)
    handle_nec(nec_geometry_complete(nec, 1))

    extype = EXTYPE
    i2 = 1 #source seg
    i3 = floor(Int64, N_SEGMENTS/2) + 1 #middle seg of odd number of segs
    i4 = 0
    tmp1 = 1.0 #excitation volts? 
    tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = 0.0
    handle_nec(nec_ex_card(nec, extype, i2, i3, i4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)) #feedpoint

    handle_nec(nec_fr_card(nec, 0, 1, FREQ_MHZ, 0.0))
    handle_nec(nec_gn_card(nec, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    handle_nec(nec_rp_card(nec, 0, 1, 1, 0, 5, 0, 0, 90.0-problem.departure_angle, 0.0, 
        0.0, 0.0, 0.0, 0.0))

    result = f(nec)
    nec_delete(nec)

    result
end

end #module
