module FlexRepresentationMethod1d

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

function assemble( degree,
                   elem_n,
                   nodes_interior,
                   q_rules_interior,
                   nodes_c_bdry,
                   nodes_t_bdry,
                   chi::Function,
                   penalty_constraint::Function,
                   constraint::Function,
                   E::Function,
                   A::Function,
                   load::Function,
                   traction::Function )

    bm_interior = BasisMesh.layout_bspline_1d_uniform_h_max_k( degree, elem_n )
    bm_interior_fc = BasisMesh.function_collection( bm_interior )
    bs_interior_fc = BasisSpline.function_collection( bm_interior_fc )
    geom_interior_fc = Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )
    mi_interior_fc = Geometry.function_collection_map_inversion_1d( bm_interior_fc, geom_interior_fc )
    q_interior = Quadrature.layout_gauss_legendre_1d( bm_interior_fc.element_count, q_rules_interior )
    q_interior_fc = Quadrature.function_collection_quadrature( q_interior )
    fs_interior_N_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc. bs_interior_fc.global_basis_value )
    fs_interior_dNdx_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_parametric_gradient )

    bm_t_bdry = BasisMesh.layout_bspline_0d()
    bm_t_bdry_fc = BasisMesh.function_collection( bm_t_bdry )
    bs_t_bdry_fc = BasisSpline.function_collection( bm_t_bdry_fc )
    geom_t_bdry_fc = Field.function_collection( bm_t_bdry_fc, bs_t_bdry_fc, nodes_t_bdry )
    mi_t_bdry_fc = Geometry.function_collection_map_inversion_1d( bm_t_bdry_fc, geom_t_bdry_fc )
    q_t_bdry = Quadrature.layout_gauss_legendre_0d()
    q_t_bdry_fc = Quadrature.function_collection_quadrature( q_t_bdry )
    fs_t_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_t_bdry_fc, mi_t_bdry_fc. bs_t_bdry_fc.global_basis_value )

    bm_c_bdry = BasisMesh.layout_bspline_0d()
    bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
    bs_c_bdry_fc = BasisSpline.function_collection( bm_c_bdry_fc )
    geom_c_bdry_fc = Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, nodes_c_bdry )
    mi_c_bdry_fc = Geometry.function_collection_map_inversion_1d( bm_c_bdry_fc, geom_c_bdry_fc )
    q_c_bdry = Quadrature.layout_gauss_legendre_0d()
    q_c_bdry_fc = Quadrature.function_collection_quadrature( q_c_bdry )
    fs_c_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc. bs_c_bdry_fc.global_basis_value )

    u_func_n = bm_interior_fc.global_function_count()
    lambda_func_n = bm_c_bdry_fc.global_function_count()

    K = zeros( u_func_n, u_func_n )
    M = zeros( u_func_n, u_func_n )
    B = zeros( lambda_func_n, u_func_n )
    F = zeros( u_func_n, 1 )
    G = zeros( lambda_func_n, 1 )
    H = zeros( u_func_n, 1 )

    weight_K( x, e, xi ) = chi( x ) * E( x ) * A( x ) * ( 1.0 / geom_interior_fc.field_parametric_gradient( e, xi ) )
    i_K_fc = Integral.function_collection_integral( q_interior_fc, geom_interior_fc, weight_K )
    K = Assemble.assembleInnerProduct!( i_K_fc, fs_interior_dNdx_fc, fs_interior_dNdx_fc, K )

    weight_M( x, e, xi ) = penalty_constraint( x )
    i_M_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_M )
    M = Assemble.assembleInnerProduct!( i_M_fc, fs_interior_N_fc, fs_interior_N_fc, M )

    weight_B( x, e, xi ) = weight_M( x, e, xi )
    i_B_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_B )
    B = Assemble.assembleInnerProduct!( i_B_fc, fs_c_bdry_N_fc, fs_interior_N_fc, B )

    weight_body_force( x, e, xi ) = load( x ) * chi( x ) * geom_interior_fc.field_parametric_gradient( e, xi )
    i_body_force_fc = Integral.function_collection_integral( q_interior_fc, geom_interior_fc, weight_body_force )
    F = Assemble.assembleProjection!( i_body_force_fc, fs_interior_N_fc, F )

    weight_traction_force( x, e, xi ) = traction( x )
    i_traction_force_fc = Integral.function_collection_integral( q_t_bdry_fc, geom_t_bdry_fc, weight_traction_force )
    F = Assemble.assembleProjection!( i_traction_force_fc, fs_interior_N_fc, F )

    weight_G( x, e, xi ) = constraint( x )
    i_G_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_G )
    G = Assemble.assembleProjection!( i_G_fc, fs_c_bdry_N_fc, G )

    weight_H( x, e, xi ) = penalty_constraint( x ) * constraint( x )
    i_H_fc = Integral.function_collection_integral( q_c_bdry_fc, geom_c_bdry_fc, weight_H )
    H = Assemble.assembleProjection!( weight_H, fs_interior_N_fc, H )

    return K, M, B, F, G, H
end

end
