import Pkg
Pkg.activate(".")
import CoreformIGA
import LinearAlgebra

#filename = "/Users/zhihui/CoreformIGA.jl/chris_old/"
#io = open(filename, "r")
#file = read(io, String)
#
#layout_interior, nodes = CoreformIGA.ImportBEXT.get1dLayoutFromBEXT( file )

function solve1dSystem( K, B, F, G )
    cN = size( B )[1]
    KB = [ [K; B] [B'; zeros( cN, cN )] ]
    Fs = [ F; G ]
    print("KB ", KB, "\n")
    print("FS ", Fs, "\n")
    #print("KB ", size(KB), "\n")
    #print("FS ", size(Fs), "\n")
    d = KB\Fs
    return d, LinearAlgebra.cond( KB )
end

filename = "/chris_old/bext_files/full_flex_bext_p2h2.json"
filepath = string(@__DIR__, filename)
io = open(filepath, "r")
file = read(io, String)
layout_interior, nodes = CoreformIGA.ImportBEXT.get1dLayoutFromBEXT( file )

nodes_interior = [ node[1] for node in nodes ]
print(nodes_interior)
quad_rules_interior( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc ) = CoreformIGA.Quadrature.layout_gauss_legendre_1d( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc, "QP1" )

geom_inv_map_interior = CoreformIGA.Geometry.function_collection_map_inversion_1d

nodes_constraint_bdry = [0]
nodes_traction_bdry = [1]

chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-9
penalty_constraint(x) = 1
constraint(x) = 0
E(x) = 1
A(x) = 1
load(x) = 0
traction(x) = 1


dirichlet_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_constraint_bdry ], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ constraint ] )
dirichlet_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( dirichlet_bc_layouts )

neumann_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_traction_bdry ], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ traction ] )
neumann_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( neumann_bc_layouts )



disp_strain_mat( x ) = CoreformIGA.Formulation1DSolid.dis_strain_mat( x )


K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod.assemble( layout_interior,
                                                                  nodes_interior,
                                                                  quad_rules_interior,
                                                                  geom_inv_map_interior,
                                                                  dirichlet_bcs_fc,
                                                                  neumann_bcs_fc,
                                                                  chi,
                                                                  penalty_constraint,
                                                                  E,
                                                                  A,
                                                                  load,
                                                                  disp_strain_mat )

println( K )
println( M )
println( B )
println( F )
println( G )
println( H )

#d, lam = CoreformIGA.NonlinearSolver.uzawaIteration( K, M, B, F, G, H )
d, cond = solve1dSystem( K, B, F, G )


println("\n", d )
