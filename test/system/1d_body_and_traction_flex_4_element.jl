using Test
using LinearAlgebra
import CoreformIGA

@testset "1d_body_and_traction_flex_4_element.jl" begin

    function solve1dSystem( K, B, F, G )
        cN = size( B )[1]
        KB = [ [K; B] [B'; zeros( cN, cN )] ]
        Fs = [ F; G ]
        d = KB\Fs
        return d, LinearAlgebra.cond( KB )
    end

    filename = "/../json/full_flex_bext_p2h2.json"
    filepath = string(@__DIR__, filename)
    io = open(filepath, "r")
    file = read(io, String)
    layout_interior, nodes = CoreformIGA.ImportBEXT.get1dLayoutFromBEXT( file )

    nodes_interior = [ node[1] for node in nodes ]
    quad_rules_interior( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc ) = CoreformIGA.Quadrature.layout_gauss_legendre_1d( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc, "QP1" )

    geom_inv_map_interior = CoreformIGA.Geometry.function_collection_map_inversion_1d

    nodes_constraint_bdry = [0]
    nodes_traction_bdry = [1]

    chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-9
    penalty_constraint(x) = 1
    constraint(x) = 0
    E(x) = 1
    A(x) = 1
    load(x) = x
    traction(x) = 1

    interior_index = zeros(Int64, 1,6)
    interior_index[1,:] = [ 1 2 3 4 5 6]
    interior_index_layout = CoreformIGA.Index.Layout( 1, 6, interior_index )
    interior_index_fc = CoreformIGA.Index.function_collection( interior_index_layout )

    bm_interior_fc = CoreformIGA.BasisMesh.function_collection( layout_interior )
    bs_interior_fc = CoreformIGA.BasisSpline.function_collection( bm_interior_fc )
    geom_interior_fc = CoreformIGA.Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )

    mi_interior_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_interior_fc, geom_interior_fc )
    q_interior_fc = CoreformIGA.Quadrature.function_collection_quadrature( quad_rules_interior( bm_interior_fc.element_count, bm_interior_fc.element_degree, mi_interior_fc, 0, 1 ) )

    fs_interior_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc )

    formulation_fc = CoreformIGA.Formulation1DSolid.function_collection()

    K_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) = chi( x ) * E( x ) * A( x ) * wi * formulation_fc.K_mat( test_fs_fc, trial_fs_fc, geom_field_fc, xi_test_i, e_test_i, xi_trial_i, e_trial_i )
    K_integral = CoreformIGA.Integral.function_collection_integral( q_interior_fc, K_integrand, "LHS_SYM" )

    K_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( K_integral,
                                                                                geom_interior_fc,
                                                                                interior_index_fc,
                                                                                interior_index_fc,
                                                                                fs_interior_fc,
                                                                                fs_interior_fc )

    F_body_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = chi( x ) * load( x ) * wi * formulation_fc.body_force_vec( test_fs_fc, geom_field_fc, xi_test_i, e_test_i )
    F_body_integral = CoreformIGA.Integral.function_collection_integral( q_interior_fc, F_body_integrand, "RHS" )
    F_body_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( F_body_integral, geom_interior_fc, interior_index_fc, interior_index_fc, fs_interior_fc, fs_interior_fc )


    # dirichlet bcs
    index = zeros(Int64,1,1)
    index[1,1] = 7
    dirichlet_index_layout = CoreformIGA.Index.Layout( 7, 7, index )
    dirichlet_index_fc = CoreformIGA.Index.function_collection( dirichlet_index_layout )

    bm_c_bdry_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_c_bdry_fc = CoreformIGA.BasisSpline.function_collection( bm_c_bdry_fc )
    geom_c_bdry_fc = CoreformIGA.Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, nodes_constraint_bdry ) 

    q_c_bdry_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_c_bdry_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_c_bdry_fc, geom_c_bdry_fc ) 
    fs_c_bdry_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc, bs_c_bdry_fc )

    B_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) =  wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i ) * ( trial_fs_fc.global_basis_evaluator( e_trial_i, xi_trial_i )' )                                    
    B_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, B_integrand, "LHS_NONSYM" )
    B_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( B_integral, geom_c_bdry_fc, dirichlet_index_fc, interior_index_fc, fs_c_bdry_fc, fs_interior_fc )

    G_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )                                    
    G_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, G_integrand, "RHS" )
    G_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( G_integral, geom_c_bdry_fc, dirichlet_index_fc, dirichlet_index_fc, fs_c_bdry_fc, fs_c_bdry_fc )

    H_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = penalty_constraint( x ) * constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )
    H_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, H_integrand, "RHS" )
    H_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( H_integral, geom_c_bdry_fc, interior_index_fc, interior_index_fc, fs_interior_fc, fs_interior_fc )

    M_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) = penalty_constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i ) * ( trial_fs_fc.global_basis_evaluator( e_trial_i, xi_trial_i )' )                                    
    M_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, M_integrand, "LHS_SYM" )
    M_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( M_integral, geom_c_bdry_fc, interior_index_fc, interior_index_fc, fs_interior_fc, fs_interior_fc )

    # Neumann bcs
    bm_t_neumann_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_t_neumann_fc = CoreformIGA.BasisSpline.function_collection( bm_t_neumann_fc )
    geom_t_neumann_fc = CoreformIGA.Field.function_collection( bm_t_neumann_fc, bs_t_neumann_fc, nodes_traction_bdry ) 

    q_t_neumann_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_t_neumann_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_t_neumann_fc, geom_t_neumann_fc ) # included in FunctionSpacce
    fs_t_neumann_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_t_neumann_fc, mi_t_neumann_fc, bs_t_neumann_fc )

    F_t_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = traction( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )
    F_t_integral = CoreformIGA.Integral.function_collection_integral( q_t_neumann_fc, F_t_integrand, "RHS" )
    F_t_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( F_t_integral, geom_t_neumann_fc, interior_index_fc, interior_index_fc, fs_interior_fc, fs_interior_fc )



    continuous_components = [ K_continuous_comp_fc, B_continuous_comp_fc, G_continuous_comp_fc, H_continuous_comp_fc, M_continuous_comp_fc, F_t_continuous_comp_fc, F_body_continuous_comp_fc ]
    K, F = CoreformIGA.Assembler.assemble( continuous_components )

    d = K\F
    cond = LinearAlgebra.cond( K )

    bm_fc = CoreformIGA.BasisMesh.function_collection( layout_interior )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    geom_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes_interior )
    field = d[1:bm_fc.global_function_count()]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, field )
    mi_interior_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_fc, geom_fc )
    x_sol = 1.0
    unused_e = 0
    unused_xi = 0
    e_inv, xi_inv = mi_interior_fc.geometric_map_inversion( x_sol, unused_e, unused_xi )
    tip_displacement = field_fc.field_value( e_inv, xi_inv )
    tol = 1e-3
    @test tip_displacement >= 4/3 - tol
    @test tip_displacement <= 4/3 + tol
end
