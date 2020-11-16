module Examples

using ..BasisSpline
using ..FlexRepresentationMethod

function stretchingAxialRod_FRM_AL( deg, elem_n, quad_rules, cad_domain, flex_domain,
                                    p_cad, p_u, E, A, load, traction, constraint )

    #deg = 2
    #elem_n = 10
    #quad_rule = 3
    #quad_rules = [ quad_rule for i in 1:elem_n ]
    #cad_domain = [ 1, 2 ]
    #flex_domain = [ 1, 3 ]
    #p_cad = 1e-12
    #p_u = 1
    #E = 1
    #A = 1
    #function load( x )
    #    return 0
    #end
    #function traction( x )
    #    return 1
    #end
    #function constraint( x )
    #    return 0
    #end
    layout = BasisSpline.buildUniformHMaxK( deg, elem_n, domain = flex_domain )
    nodes = BasisSpline.nodesEquallySpaced( layout )
    return FlexRepresentationMethod1d.solve( layout, nodes, quad_rules,
                                             p_cad, cad_domain, E, A, load,
                                             traction, p_u, constraint )
end

end
