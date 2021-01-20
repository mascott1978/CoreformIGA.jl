using Test
using LinearAlgebra
import CoreformIGA

@testset "1d_traction_1element.jl" begin
    layout_interior = CoreformIGA.CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 1 )
    nodes_interior = [ [0, 0, 0, 1] [1, 0, 0, 1] ]
    quad_rules_interior( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc ) = CoreformIGA.Quadrature.layout_gauss_legendre_1d( elem_count, elem_deg, inverse_map, d_bc_fc, t_bc_fc, "QP1" )
    geom_inv_map_interior = 
    nodes_constraint_bdry = zeros(4,1)
    nodes_constraint_bdry[ 4 ] = 1
    nodes_traction_bdry = zeros(4,1)
    nodes_traction_bdry[ 1 ] = 1
    nodes_traction_bdry[ 4 ] = 1
    chi(x) = x[1] >= 0 && x[1] <= 1 ? 1.0 : 1e-12
    penalty_constraint(x) = 1
    constraint(x) = 0
    E(x) = 1
    A(x) = 1
    load(x) = 0
    traction(x) = 1

    # stiffness
    interior_index = zeros(Int64, 1,2)
    interior_index[1,:] = [ 1; 2 ]
    interior_index_layout = CoreformIGA.Index.Layout( 1, 2, interior_index )
    interior_index_fc = CoreformIGA.Index.function_collection( interior_index_layout )

    bm_interior_fc = CoreformIGA.BasisMesh.function_collection( layout_interior )
    bs_interior_fc = CoreformIGA.BasisSpline.function_collection( bm_interior_fc )
    geom_interior_fc = CoreformIGA.Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )

    mi_interior_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_interior_fc, geom_interior_fc )
    q_interior_fc = CoreformIGA.Quadrature.function_collection_quadrature( quad_rules_interior( bm_interior_fc.element_count, bm_interior_fc.element_degree, mi_interior_fc, 0, 1 ) )

    fs_interior_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc )

    formulation_fc = CoreformIGA.Formulation1DSolid.function_collection()

    K_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) = chi( x ) * E( x ) * A( x ) * wi * formulation_fc.K_mat( test_fs_fc, trial_fs_fc, geom_field_fc, xi_test_i, e_test_i, xi_trial_i, e_trial_i )
    K_integral = CoreformIGA.Integral.function_collection_integral( q_interior_fc, K_integrand, CoreformIGA.Integral.LHS_SYM )
    # K_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( K_integral,
    #                                                                             geom_interior_fc,
    #                                                                             interior_index_fc,
    #                                                                             interior_index_fc,
    #                                                                             fs_interior_fc,
    #                                                                             fs_interior_fc )
    K_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral = K_integral,
                                                                                geom_field = geom_interior_fc,
                                                                                test_index = interior_index_fc,
                                                                                trial_index = interior_index_fc,
                                                                                test_function_space = fs_interior_fc,
                                                                                trial_function_space = fs_interior_fc )

    # body force
    F_body_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = chi( x ) * load( x ) * wi * formulation_fc.body_force_vec( test_fs_fc, geom_field_fc, xi_test_i, e_test_i )
    F_body_integral = CoreformIGA.Integral.function_collection_integral( q_interior_fc, F_body_integrand, CoreformIGA.Integral.RHS )
    F_body_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = F_body_integral,
                                                                                     geom_field           = geom_interior_fc,
                                                                                     test_index           = interior_index_fc,
                                                                                     trial_index          = interior_index_fc,
                                                                                     test_function_space  = fs_interior_fc,
                                                                                     trial_function_space = fs_interior_fc )


    # dirichlet bcs
    index = zeros(Int64,1,1)
    index[1,1] = 3
    dirichlet_index_layout = CoreformIGA.Index.Layout( 3, 3, index )
    dirichlet_index_fc = CoreformIGA.Index.function_collection( dirichlet_index_layout )

    bm_c_bdry_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_c_bdry_fc = CoreformIGA.BasisSpline.function_collection( bm_c_bdry_fc )
    geom_c_bdry_fc = CoreformIGA.Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, nodes_constraint_bdry ) 

    q_c_bdry_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_c_bdry_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_c_bdry_fc, geom_c_bdry_fc ) 
    fs_c_bdry_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc, bs_c_bdry_fc )

    B_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) =  wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i ) * ( trial_fs_fc.global_basis_evaluator( e_trial_i, xi_trial_i )' )                                    
    B_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, B_integrand, CoreformIGA.Integral.LHS_NONSYM )
    B_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = B_integral,
                                                                                geom_field           = geom_c_bdry_fc,
                                                                                test_index           = dirichlet_index_fc,
                                                                                trial_index          = interior_index_fc,
                                                                                test_function_space  = fs_c_bdry_fc,
                                                                                trial_function_space = fs_interior_fc )

    G_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )                                    
    G_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, G_integrand, CoreformIGA.Integral.RHS )
    G_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = G_integral,
                                                                                geom_field           = geom_c_bdry_fc,
                                                                                test_index           = dirichlet_index_fc,
                                                                                trial_index          = dirichlet_index_fc,
                                                                                test_function_space  = fs_c_bdry_fc,
                                                                                trial_function_space = fs_c_bdry_fc )

    H_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = penalty_constraint( x ) * constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )
    H_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, H_integrand, CoreformIGA.Integral.RHS )
    H_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = H_integral,
                                                                                geom_field           = geom_c_bdry_fc,
                                                                                test_index           = interior_index_fc,
                                                                                trial_index          = interior_index_fc,
                                                                                test_function_space  = fs_interior_fc,
                                                                                trial_function_space = fs_interior_fc )

    M_integrand( x, xi_test_i, e_test_i, xi_trial_i, e_trial_i, wi, geom_field_fc, test_fs_fc, trial_fs_fc ) = penalty_constraint( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i ) * ( trial_fs_fc.global_basis_evaluator( e_trial_i, xi_trial_i )' )                                    
    M_integral = CoreformIGA.Integral.function_collection_integral( q_c_bdry_fc, M_integrand, CoreformIGA.Integral.LHS_SYM )
    M_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = M_integral,
                                                                                geom_field           = geom_c_bdry_fc,
                                                                                test_index           = interior_index_fc,
                                                                                trial_index          = interior_index_fc,
                                                                                test_function_space  = fs_interior_fc,
                                                                                trial_function_space = fs_interior_fc )

    # Neumann bcs
    bm_t_neumann_fc = CoreformIGA.BasisMesh.function_collection( CoreformIGA.BasisMesh.layout_bspline_0d() )
    bs_t_neumann_fc = CoreformIGA.BasisSpline.function_collection( bm_t_neumann_fc )
    geom_t_neumann_fc = CoreformIGA.Field.function_collection( bm_t_neumann_fc, bs_t_neumann_fc, nodes_traction_bdry ) 

    q_t_neumann_fc = CoreformIGA.Quadrature.function_collection_quadrature( CoreformIGA.Quadrature.layout_gauss_legendre_0d() )

    mi_t_neumann_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_t_neumann_fc, geom_t_neumann_fc ) # included in FunctionSpacce
    fs_t_neumann_fc = CoreformIGA.FunctionSpace.function_collection_function_space( bm_t_neumann_fc, mi_t_neumann_fc, bs_t_neumann_fc )

    F_t_integrand( x, xi_test_i, e_test_i, wi, geom_field_fc, test_fs_fc ) = traction( x ) * wi * test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )
    F_t_integral = CoreformIGA.Integral.function_collection_integral( q_t_neumann_fc, F_t_integrand, CoreformIGA.Integral.RHS )
    F_t_continuous_comp_fc = CoreformIGA.ContinuousComponent.function_collection( integral             = F_t_integral,
                                                                                  geom_field           = geom_t_neumann_fc,
                                                                                  test_index           = interior_index_fc,
                                                                                  trial_index          = interior_index_fc,
                                                                                  test_function_space  = fs_interior_fc,
                                                                                  trial_function_space = fs_interior_fc )



    # neumann_bc_fc = CoreformIGA.ContinuousComponent.function_collection( traction, dirichlet_index_fc, geom_t_neumann_fc, q_t_neumann_fc, fs_t_neumann_fc )

    continuous_components = [ K_continuous_comp_fc, B_continuous_comp_fc, G_continuous_comp_fc, H_continuous_comp_fc, M_continuous_comp_fc, F_t_continuous_comp_fc, F_body_continuous_comp_fc ]
    K, F = CoreformIGA.Assembler.assemble( continuous_components )

    # F_ref = zeros( 2, 1 );
    # F_ref[2] = 1.0;
    # G_ref = zeros( 1, 1 );
    # H_ref = zeros( 2, 1 );
    # @test K == [1.0 -1.0; -1.0 1.0]
    # @test M == [1.0 0.0; 0.0 0.0]
    # @test B == [1.0 0.0]
    # @test F == F_ref
    # @test G == G_ref
    # @test H == H_ref

    F_ref = zeros( 3, 1 )
    F_ref[ 2 ] = 1
    K_ref = ones( 3, 3 )
    K_ref[ 1, 1 ] = 2
    K_ref[ 1, 2 ] = -1
    K_ref[ 1, 2 ] = -1
    K_ref[ 2, 1 ] = -1
    K_ref[ 2, 3 ] = 0
    K_ref[ 3, 2 ] = 0
    K_ref[ 3, 3 ] = 0
    @test F == F_ref
    @test K == K_ref
end
