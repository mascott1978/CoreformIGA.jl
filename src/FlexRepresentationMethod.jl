module FlexRepresentationMethod

using ..NonlinearSolver
using ..Assembler
using ..Quadrature
using ..BasisBernstein
using ..BasisMesh
using ..BasisSpline
using ..FunctionSpace
using ..Field
using ..Geometry
using ..Integral

function assemble( bm_interior, # Why are we calling this interior? This seems like the flex mesh
                   nodes_interior,
                   q_rules_interior,
                   geom_inv_map_interior,
                   dirichelet_bcs,
                   neumann_bcs,
                   chi::Function,
                   penalty_constraint::Function,
                   E::Function,
                   A::Function,
                   load::Function,
                   disp_strain_mat::Function
                   )

    #bm_interior = BasisMesh.layout_bspline_1d_uniform_h_max_k( degree, elem_n )
    bm_interior_fc = BasisMesh.function_collection( bm_interior )
    bs_interior_fc = BasisSpline.function_collection( bm_interior_fc )
    geom_interior_fc = Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )
    mi_interior_fc = geom_inv_map_interior( bm_interior_fc, geom_interior_fc )
    q_interior = q_rules_interior( bm_interior_fc.element_count, bm_interior_fc.element_degree, mi_interior_fc, dirichelet_bcs, neumann_bcs )
    q_interior_fc = Quadrature.function_collection_quadrature( q_interior )
    fs_interior_N_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_value )
    fs_interior_dNdx_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_parametric_gradient )
    u_func_n = bm_interior_fc.global_function_count()
    lambda_func_n = 0
    for i = 1:dirichelet_bcs.bc_num()
        bm_c_bdry = dirichelet_bcs.basis_mesh( i )
        bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
        lambda_func_n = lambda_func_n + bm_c_bdry_fc.global_function_count()
    end
    K = zeros( u_func_n, u_func_n )
    M = zeros( u_func_n, u_func_n )
    B = zeros( lambda_func_n, u_func_n )
    F = zeros( u_func_n, 1 )
    G = zeros( lambda_func_n, 1 )
    H = zeros( u_func_n, 1 )

    # Dirichlet BCs
    for i = 1:dirichelet_bcs.bc_num()
        bm_c_bdry = dirichelet_bcs.basis_mesh( i )
        bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
        bs_c_bdry_fc = BasisSpline.function_collection( bm_c_bdry_fc )
        geom_c_bdry_fc = Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, dirichelet_bcs.nodes( i ) )
        mi_c_bdry_fc = dirichelet_bcs.geom_inv_map( i )( bm_c_bdry_fc, geom_c_bdry_fc )
        q_c_bdry = dirichelet_bcs.quadrature_layout( i )
        q_c_bdry_fc = Quadrature.function_collection_quadrature( q_c_bdry )
        fs_c_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc, bs_c_bdry_fc.global_basis_value )
        lambda_func_n = lambda_func_n + bm_c_bdry_fc.global_function_count()
        weight_B( x, e, xi ) = 1
        i_B_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_B )
        B = Assembler.assembleInnerProduct!( i_B_fc, fs_c_bdry_N_fc, fs_interior_N_fc, disp_strain_mat, B ) # FIXME: disp_strain_mat might be a different function for 2D

        weight_G( x, e, xi ) = dirichelet_bcs.value( i )( x )
        i_G_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_G )
        G = Assembler.assembleProjection!( i_G_fc, fs_c_bdry_N_fc, G )

        weight_H( x, e, xi ) = penalty_constraint( x ) * dirichelet_bcs.value( i )( x )
        i_H_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_H )
        H = Assembler.assembleProjection!( i_H_fc, fs_interior_N_fc, H )

        weight_M( x, e, xi ) = penalty_constraint( x )
        i_M_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_M )
        M = Assembler.assembleInnerProduct!( i_M_fc, fs_interior_N_fc, fs_interior_N_fc, disp_strain_mat, M ) # FIXME: disp_strain_mat might be a different function for 2D

    end

    # Neumann Bcs
    for i = 1:neumann_bcs.bc_num()
        bm_t_bdry = neumann_bcs.basis_mesh( i )
        bm_t_bdry_fc = BasisMesh.function_collection( bm_t_bdry )
        bs_t_bdry_fc = BasisSpline.function_collection( bm_t_bdry_fc )
        geom_t_bdry_fc = Field.function_collection( bm_t_bdry_fc, bs_t_bdry_fc, neumann_bcs.nodes( i ) )
        mi_t_bdry_fc = neumann_bcs.geom_inv_map( i )( bm_t_bdry_fc, geom_t_bdry_fc )
        q_t_bdry = neumann_bcs.quadrature_layout( i )
        q_t_bdry_fc = Quadrature.function_collection_quadrature( q_t_bdry )
        fs_t_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_t_bdry_fc, mi_t_bdry_fc, bs_t_bdry_fc.global_basis_value )

        weight_traction_force( x, e, xi ) = neumann_bcs.value( i )( x )
        i_traction_force_fc = Integral.function_collection_integral( q_t_bdry_fc, geom_t_bdry_fc, weight_traction_force )
        F = Assembler.assembleProjection!( i_traction_force_fc, fs_interior_N_fc, F )
    end

    weight_K( x, e, xi ) = chi( x ) * E( x ) * A( x ) * ( 1.0 / geom_interior_fc.field_parametric_gradient( e, xi ) )
    i_K_fc = Integral.function_collection_integral( q_interior_fc, geom_interior_fc, weight_K )
    K = Assembler.assembleInnerProduct!( i_K_fc, fs_interior_dNdx_fc, fs_interior_dNdx_fc, disp_strain_mat, K )

    weight_body_force( x, e, xi ) = load( x ) * chi( x ) * geom_interior_fc.field_parametric_gradient( e, xi )
    i_body_force_fc = Integral.function_collection_integral( q_interior_fc, geom_interior_fc, weight_body_force )
    F = Assembler.assembleProjection!( i_body_force_fc, fs_interior_N_fc, F )

    return K, M, B, F, G, H
end

end
