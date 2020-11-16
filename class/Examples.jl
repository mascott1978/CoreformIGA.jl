
# stretching of an axial rod, using FRM AL

deg = 2
elem_n = 20
quad_rule = 7
quad_rules = [ quad_rule for i in 1:elem_n ]
cad_domain = [ 1, 2 ]
flex_domain = [ 1, 4 ]
p_cad = 1e-12
p_u = 1
E = 1
A = 1
function load( x )
    return 0
end
function traction( x )
    return 1
end
function constraint( x )
    return 0
end
layout = CoreformIGA.BasisSpline.buildUniformHMaxK( deg, elem_n, domain = flex_domain )
nodes = CoreformIGA.BasisSpline.nodesEquallySpaced( layout )
d, lambda = CoreformIGA.FlexRepresentationMethod1d.solve( layout, nodes, quad_rules,
                                                          p_cad, cad_domain, E, A, load,
                                                          traction, p_u, constraint )
import Plots
plt = Plots.plot()
for e in 1:elem_n
    CoreformIGA.Field.plot!( plt, layout, nodes, d, e, CoreformIGA.Field.F )
end
plt
