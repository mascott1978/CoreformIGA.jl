using Test
using LinearAlgebra
import CoreformIGA

@testset "1d_traction_1element.jl" begin
    layout_interior = CoreformIGA.CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 1 )
    nodes_interior = [0,1]
    quad_rules_interior( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc ) = CoreformIGA.Quadrature.layout_gauss_legendre_1d( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc, "QP1" )
    geom_inv_map_interior = 
    nodes_constraint_bdry = [0]
    nodes_traction_bdry = [1]
    chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12
    penalty_constraint(x) = ones(1,1)
    constraint(x) = zeros(1,1)
    E(x) = 1
    A(x) = 1
    load(x) = 0
    traction(x) = ones( 1, 1 )

    # dirichlet_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_constraint_bdry ], [ [ [ 1 ] ] ],  [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ constraint ] )
    # dirichlet_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( dirichlet_bc_layouts )

    index = zeros(Int64,1,1)
    index[1,1] = 1
    dirichlet_index_layout = CoreformIGA.Index.Layout( 1, 1, index )
    dirichlet_index_fc = CoreformIGA.Index.function_collection( dirichlet_index_layout )

    bm_c_bdry_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_c_bdry_fc = CoreformIGA.BasisSpline.function_collection( bm_c_bdry_fc )
    geom_c_bdry_fc = CoreformIGA.Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, nodes_constraint_bdry ) 

    q_c_bdry_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_c_bdry_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_c_bdry_fc, geom_c_bdry_fc ) 
    fs_c_bdry_N_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc, bs_c_bdry_fc.global_basis_value )

    dirichlet_bc_fc = CoreformIGA.ContinuousComponent.function_collection( constraint, dirichlet_index_fc, geom_c_bdry_fc, q_c_bdry_fc, fs_c_bdry_N_fc )

    # neumann_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [ nodes_traction_bdry ], [ [ [ 1 ] ] ], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ traction ] )
    # neumann_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( neumann_bc_layouts )


    bm_t_neumann_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_t_neumann_fc = CoreformIGA.BasisSpline.function_collection( bm_t_neumann_fc )
    geom_t_neumann_fc = CoreformIGA.Field.function_collection( bm_t_neumann_fc, bs_t_neumann_fc, nodes_traction_bdry ) 

    q_t_neumann_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_t_neumann_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_t_neumann_fc, geom_t_neumann_fc ) # included in FunctionSpacce
    fs_t_neumann_N_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_t_neumann_fc, mi_t_neumann_fc, bs_t_neumann_fc.global_basis_value )

    neumann_bc_fc = CoreformIGA.ContinuousComponent.function_collection( traction, dirichlet_index_fc, geom_t_neumann_fc, q_t_neumann_fc, fs_t_neumann_N_fc )


    interior_index = zeros(Int64, 1,2)
    interior_index[1,:] = [ 1 2]
    interior_index_layout = CoreformIGA.Index.Layout( 1, 2, interior_index )
    interior_index_fc = CoreformIGA.Index.function_collection( interior_index_layout )

    bm_interior_fc = CoreformIGA.BasisMesh.function_collection( layout_interior )
    bs_interior_fc = CoreformIGA.BasisSpline.function_collection( bm_interior_fc )
    geom_interior_fc = CoreformIGA.Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )

    mi_interior_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_interior_fc, geom_interior_fc )
    q_interior_fc = CoreformIGA.Quadrature.function_collection_quadrature( quad_rules_interior( bm_interior_fc.element_count, bm_interior_fc.element_degree, mi_interior_fc, dirichlet_bc_fc, neumann_bc_fc ) )

    fs_interior_N_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_value )
    fs_interior_dNdx_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_parametric_gradient )

    interior_fc = CoreformIGA.ContinuousComponent.function_collection( load, interior_index_fc, geom_interior_fc, q_interior_fc, fs_interior_N_fc, fs_interior_dNdx_fc )

    formulation = CoreformIGA.Formulation1DSolid.function_collection()

    K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod.assemble( interior_fc,
                                                                      [ dirichlet_bc_fc ],
                                                                      [ neumann_bc_fc ],
                                                                      formulation,
                                                                      chi,
                                                                      penalty_constraint,
                                                                      E,
                                                                      A )
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
