
# stretching of an axial rod, using FRM AL


import CoreformIGA
layout_interior = CoreformIGA.BasisMesh.bspline( 1, 1 )
nodes_interior = [0,1]
quad_rules_interior = [2]

layout_constraint_bdry = CoreformIGA.BasisMesh.vertex()
nodes_constraint_bdry = [0]

layout_traction_bdry = CoreformIGA.BasisMesh.vertex()
nodes_traction_bdry = [1]

chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12
penalty_constraint(x) = 1
constraint(x) = 0
E(x) = 1
A(x) = 1
load(x) = 0
traction(x) = 1

K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod1d.assemble( layout_interior, nodes_interior, quad_rules_interior,
                                                                    layout_constraint_bdry, nodes_constraint_bdry,
                                                                    layout_traction_bdry, nodes_traction_bdry,
                                                                    chi, penalty_constraint, constraint, E, A, load, traction)

