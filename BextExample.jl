
import CoreformIGA

filename = "/home/chris/cf/master/iga/chris_old/bext_files/full_flex_bext_p2h2.json"
io = open(filename, "r")
file = read(io, String)

layout_interior, nodes = CoreformIGA.ImportBEXT.get1dLayoutFromBEXT( file )

nodes_interior = [ node[1] for node in nodes ]
quad_rules_interior = [] #NOTE UNUSED


#layout_constraint_bdry = CoreformIGA.BasisMesh.layout_bspline_0d()
nodes_constraint_bdry = [0]

#layout_traction_bdry = CoreformIGA.BasisMesh.layout_bspline_0d()
nodes_traction_bdry = [1]

chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12





#layout_interior = CoreformIGA.BasisMesh.bspline( 1, 1 )
#nodes_interior = [0,1]
#quad_rules_interior = [2]
#
#
#layout_traction_bdry = CoreformIGA.BasisMesh.vertex()
#nodes_traction_bdry = [1]
#
penalty_constraint(x) = 1
constraint(x) = 0
E(x) = 1
A(x) = 1
load(x) = 0
traction(x) = 1
#
K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod1d.assemble( layout_interior, nodes_interior, quad_rules_interior,
                                                                    nodes_constraint_bdry,
                                                                    nodes_traction_bdry,
                                                                    chi, penalty_constraint, constraint, E, A, load, traction)
