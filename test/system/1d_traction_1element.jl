using Test
using LinearAlgebra
import CoreformIGA

@testset "1d_traction_1element.jl" begin
    layout_interior = CoreformIGA.CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 1 )
    nodes_interior = [0,1]
    quad_rules_interior( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc ) = CoreformIGA.Quadrature.layout_gauss_legendre_1d( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc, "QP1" )
    geom_inv_map_interior = CoreformIGA.Geometry.function_collection_map_inversion_1d
    nodes_constraint_bdry = [0]
    nodes_traction_bdry = [1]
    chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12
    penalty_constraint(x) = ones(1,1)
    constraint(x) = zeros(1,1)
    E(x) = 1
    A(x) = 1
    load(x) = 0
    traction(x) = ones( 1, 1 )

    dirichlet_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_constraint_bdry ], [ [ [ 1 ] ] ],  [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ constraint ] )
    dirichlet_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( dirichlet_bc_layouts )

    neumann_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_traction_bdry ], [ [ [ 1 ] ] ], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ traction ] )
    neumann_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( neumann_bc_layouts )

    formulation = CoreformIGA.Formulation1DSolid.function_collection()
    K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod.assemble( layout_interior,
                                                                      nodes_interior,
                                                                      quad_rules_interior,
                                                                      geom_inv_map_interior,
                                                                      dirichlet_bcs_fc,
                                                                      neumann_bcs_fc,
                                                                      formulation,
                                                                      chi,
                                                                      penalty_constraint,
                                                                      E,
                                                                      A,
                                                                      load )

    F_ref = zeros( 2, 1 );
    F_ref[2] = 1.0;
    G_ref = zeros( 1, 1 );
    H_ref = zeros( 2, 1 );
    @test K == [1.0 -1.0; -1.0 1.0]
    @test M == [1.0 0.0; 0.0 0.0]
    @test B == [1.0 0.0]
    @test F == F_ref
    @test G == G_ref
    @test H == H_ref

end
