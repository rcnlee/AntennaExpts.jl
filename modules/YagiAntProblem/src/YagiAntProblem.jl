"""
Yagi-like antenna fitting inside a Bott's dot
"""
module YagiAntProblem

export YagiAnt, create_grammar, get_grammar, get_fitness

using ExprSearch
import ExprSearch: ExprProblem, get_fitness, get_grammar
using RLESUtils, Interpreter
using NECPP

#all measurements in meters
const MAXLENGTH = 1.0
const MAXWIDTH = 1.0
const HEIGHT = 0.2
const N_SEGMENTS = 11
const WIRE_RADIUS = 0.005 #meters
const EXTYPE = 0 #0=voltage source
const EX_TAGID = 1 #source tagid
const TAGID = 0 #default tagid

const DIVCONST = 10.0

type YagiAnt <: ExprProblem
    departure_angle::Float64 #degrees rel horizontal
    grammar::Grammar
end

function YagiAnt(departure_angle::Float64=10.0)
    grammar = create_grammar()
    YagiAnt(departure_angle, grammar)
end

function create_grammar()
    @grammar grammar begin
        start = ant
        ant[make_block] = base_el + exc_el + addl_el
        base_el = Expr(:call, :base_element, :nec, frac)
        exc_el = straight_exc_el | bent_exc_el
        straight_exc_el = Expr(:call, :straight_element, :nec, frac, frac, true) 
        bent_exc_el = Expr(:call, :bent_element, :nec, frac, frac, frac, dir, true)
        addl_el = (el)^(1:5) 
        el = straight_el | bent_el
        straight_el = Expr(:call, :straight_element, :nec, frac, frac, false)
        bent_el = Expr(:call, :bent_element, :nec, frac, frac, frac, dir, false)
        frac[divconst] = 0:10 #fraction of some length
        dir = 1.0 | -1.0
    end
    grammar
end

###
#functions for grammar
make_block(lst::Array) = Expr(:block, lst...)

function base_element(nec::NecContext, base_len_frac::Float64)
    len = base_len_frac * MAXLENGTH
    tag_id = TAG_ID 
    segment_count = N_SEGMENTS
    xw1 = 0.0
    yw1 = 0.0
    zw1 = HEIGHT 
    xw2 = 0.0
    yw2 = len 
    zw2 = HEIGHT
    rad = WIRE_RADIUS
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad))
end

function straight_element(nec::NecContext, offset_frac::Float64, width_frac::Float64, b_exc::Bool)
    offset = offset_frac * MAXLENGTH
    halfwidth = width_frac * MAXWIDTH / 2
    tag_id = b_exc ? EXC_TAGID : TAG_ID
    segment_count = N_SEGMENTS
    xw1 = -halfwidth 
    yw1 = offset 
    zw1 = HEIGHT 
    xw2 = +halfwidth 
    yw2 = offset 
    zw2 = HEIGHT
    rad = WIRE_RADIUS
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad))
end

function bent_element(nec::NecContext, offset_frac::Float64, width_frac::Float64,
    len_frac::Float64, dir::Float64, b_exc::Bool)
    offset = offset_frac * MAXLENGTH
    halfwidth = width_frac * MAXWIDTH / 2
    len = len_frac * MAXLENGTH

    tag_id = b_exc ? EX_TAGID : TAG_ID
    segment_count = N_SEGMENTS
    xw1 = -halfwidth 
    yw1 = offset 
    zw1 = HEIGHT 
    xw2 = +halfwidth 
    yw2 = offset 
    zw2 = HEIGHT
    rad = WIRE_RADIUS
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad)) #base

    tag_id = TAG_ID 
    segment_count = N_SEGMENTS
    xw1 = -halfwidth 
    yw1 = offset
    zw1 = HEIGHT 
    xw2 = -halfwidth 
    yw2 = min(MAXLENGTH, max(0.0, offset+dir*len)) #cap, don't exceed boundary
    zw2 = HEIGHT
    rad = WIRE_RADIUS
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad)) #left bend

    tag_id = TAG_ID 
    segment_count = N_SEGMENTS
    xw1 = +halfwidth 
    yw1 = offset 
    zw1 = HEIGHT 
    xw2 = +halfwidth 
    yw2 = min(MAXLENGTH, max(0.0, offset+dir*len)) #cap, don't exceed boundary
    zw2 = HEIGHT
    rad = WIRE_RADIUS
    rdel = 1.0 #no tapering
    rrad = 1.0 #no tapering
    handle_nec(nec_wire(nec, tag_id, segment_count, xw1, yw1, zw1,
        xw2, yw2, zw2, rad, rdel, rrad)) #right bend
end

divconst(x) = x / DIVCONST
###

const SYMTABLE = SymbolTable(
    :make_block => make_block,
    :base_element => base_element,
    :straight_element => straight_element,
    :bent_element => bent_element,
    :divconst => divconst
    )

function eval_expr(problem::YagiAnt, expr, nec) 
    SYMTABLE[:nec] = nec
    interpret(SYMTABLE, expr)
end

ExprSearch.get_grammar(problem::YagiAnt) = problem.grammar

function ExprSearch.get_fitness(problem::YagiAnt, expr)

    #call nec on expr here
    nec = nec_create()
    eval_expr(problem, expr, nec)
    handle_nec(nec_geometry_complete(nec, 1))

    extype = EXTYPE
    i2 = EX_TAGID #source seg
    i3 = floor(Int64, N_SEGMENTS/2) + 1 #middle seg of odd number of segs
    tmp1 = 1.0 #exc volts? 
    tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = 0.0
    handle_nec(nec_ex_card(nec, extype, i2, i3, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)) #feedpoint

    handle_nec(nec_gn_card(nec, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    handle_nec(nec_rp_card(nec, 0, 90, 1, 0, 5, 0, 0, 0.0, 90.0, 1.0, 0.0, 0.0, 0.0))

    #for now... look at impedance
    result_index = 0
    z = Complex(nec_impedance_real(nec,result_index), nec_impedance_imag(nec,result_index))
    fitness = abs(z) #magnitude 

    nec_delete(nec)

    fitness
end

end #module
