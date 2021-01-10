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
using ..Formulation1DSolid

function assemble( interior_fc,
                   # bm_interior, # Why are we calling this interior? This seems like the flex mesh
                   # nodes_interior,
                   # q_rules_interior,
                   # geom_inv_map_interior,
                   dirichlet_bcs,
                   neumann_bcs,
                   formulation,
                   chi::Function,
                   penalty_constraint::Function,
                   E::Function,
                   A::Function )

    # #bm_interior = BasisMesh.layout_bspline_1d_uniform_h_max_k( degree, elem_n )
    # bm_interior_fc = BasisMesh.function_collection( bm_interior )
    # bs_interior_fc = BasisSpline.function_collection( bm_interior_fc )
    # geom_interior_fc = Field.function_collection( bm_interior_fc, bs_interior_fc, nodes_interior )
    # mi_interior_fc = geom_inv_map_interior( bm_interior_fc, geom_interior_fc )
    # q_interior = q_rules_interior( bm_interior_fc.element_count, bm_interior_fc.element_degree, mi_interior_fc, dirichlet_bcs, neumann_bcs )
    # q_interior_fc = Quadrature.function_collection_quadrature( q_interior )
    # fs_interior_N_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_value )
    # fs_interior_dNdx_fc = FunctionSpace.function_collection_function_space( bm_interior_fc, mi_interior_fc, bs_interior_fc.global_basis_parametric_gradient )

    # # pass basis information into functions to build this.
    # u_func_n = bm_interior_fc.global_function_count()
    # node_dofn = formulation.dofn()
    # u_func_id2dof_ids = zeros( Int, u_func_n, node_dofn )
    # u_dof_n = 0
    # for i = 1:u_func_n
    #     u_func_id2dof_ids[ i, : ] = [ 1:node_dofn; ] .+ u_dof_n
    #     u_dof_n = u_dof_n + node_dofn
    # end

    # lambda_func_n = 0
    # lambda_dof_n = 0
    # lambda_func_id2dof_ids = Array{Any}( undef, dirichlet_bcs.bc_num() )
    # for i = 1:dirichlet_bcs.bc_num()
    #     bm_c_bdry = dirichlet_bcs.basis_mesh( i )
    #     bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
    #     i_lambda_func_n = bm_c_bdry_fc.global_function_count()
    #     i_dofs = dirichlet_bcs.dofs( i )
    #     lambda_func_id2dof_ids[ i ] = zeros( Int, i_lambda_func_n, node_dofn )
    #     for j = 1:i_lambda_func_n
    #         node_i_constraint_num = length( i_dofs[ j ] )
    #         lambda_func_id2dof_ids[ i ][ j, i_dofs[ j ] ] = [ 1:node_i_constraint_num; ] .+ lambda_dof_n
    #         lambda_dof_n = lambda_dof_n + node_i_constraint_num
    #     end
    #     lambda_func_n = lambda_func_n + i_lambda_func_n
    # end
    u_dof_n = interior_fc.index_fc.end_id()
    lambda_dof_n = dirichlet_bcs[ end ].index_fc.end_id()
    K = zeros( u_dof_n, u_dof_n )
    M = zeros( u_dof_n, u_dof_n )
    B = zeros( lambda_dof_n, u_dof_n )
    F = zeros( u_dof_n, 1 )
    G = zeros( lambda_dof_n, 1 )
    H = zeros( u_dof_n, 1 )

    # Dirichlet BCs
    for i = 1:length( dirichlet_bcs )
        # bm_c_bdry = dirichlet_bcs.basis_mesh( i )
        # bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
        # bs_c_bdry_fc = BasisSpline.function_collection( bm_c_bdry_fc )
        # geom_c_bdry_fc = Field.function_collection( bm_c_bdry_fc, bs_c_bdry_fc, dirichlet_bcs.nodes( i ) ) # included in Field
        # mi_c_bdry_fc = dirichlet_bcs.geom_inv_map( i )( bm_c_bdry_fc, geom_c_bdry_fc ) # included in FunctionSpacce
        # q_c_bdry = dirichlet_bcs.quadrature_layout( i )
        # q_c_bdry_fc = Quadrature.function_collection_quadrature( q_c_bdry )
        # fs_c_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_c_bdry_fc, mi_c_bdry_fc, bs_c_bdry_fc.global_basis_value )
        weight_B( x, e, xi ) = 1
        i_B_fc = Integral.function_collection_integral( dirichlet_bcs[ i ].quadrature_fc, dirichlet_bcs[ i ].geom_field_fc, weight_B )
        B = Assembler.assembleInnerProduct!( i_B_fc, dirichlet_bcs[ i ].fs_N_fc, interior_fc.fs_N_fc, formulation.dis_strain_mat, dirichlet_bcs[ i ].index_fc, interior_fc.index_fc, B ) 

        weight_G( x, e, xi ) = dirichlet_bcs[ i ].f( x )
        i_G_fc = Integral.function_collection_integral( dirichlet_bcs[ i ].quadrature_fc, dirichlet_bcs[ i ].geom_field_fc, weight_G )
        G = Assembler.assembleProjection!( i_G_fc, dirichlet_bcs[ i ].fs_N_fc, dirichlet_bcs[ i ].index_fc, G )

        weight_H( x, e, xi ) = penalty_constraint( x ) * dirichlet_bcs[ i ].f( x )
        i_H_fc = Integral.function_collection_integral( dirichlet_bcs[ i ].quadrature_fc, dirichlet_bcs[ i ].geom_field_fc, weight_H )
        H = Assembler.assembleProjection!( i_H_fc, interior_fc.fs_N_fc, interior_fc.index_fc, H )

        weight_M( x, e, xi ) = penalty_constraint( x )
        i_M_fc = Integral.function_collection_integral( dirichlet_bcs[ i ].quadrature_fc, dirichlet_bcs[ i ].geom_field_fc, weight_M )
        M = Assembler.assembleInnerProduct!( i_M_fc, interior_fc.fs_N_fc, interior_fc.fs_N_fc, formulation.dis_strain_mat, interior_fc.index_fc, interior_fc.index_fc, M )
    end

    # Neumann Bcs
    for i = 1:length( neumann_bcs )
        # bm_t_bdry = neumann_bcs.basis_mesh( i )
        # bm_t_bdry_fc = BasisMesh.function_collection( bm_t_bdry )
        # bs_t_bdry_fc = BasisSpline.function_collection( bm_t_bdry_fc )
        # geom_t_bdry_fc = Field.function_collection( bm_t_bdry_fc, bs_t_bdry_fc, neumann_bcs.nodes( i ) )
        # mi_t_bdry_fc = neumann_bcs.geom_inv_map( i )( bm_t_bdry_fc, geom_t_bdry_fc )
        # q_t_bdry = neumann_bcs.quadrature_layout( i )
        # q_t_bdry_fc = Quadrature.function_collection_quadrature( q_t_bdry )
        # fs_t_bdry_N_fc = FunctionSpace.function_collection_function_space( bm_t_bdry_fc, mi_t_bdry_fc, bs_t_bdry_fc.global_basis_value )

        weight_traction_force( x, e, xi ) = neumann_bcs[ i ].f( x )
        i_traction_force_fc = Integral.function_collection_integral( neumann_bcs[ i ].quadrature_fc, neumann_bcs[ i ].geom_field_fc, weight_traction_force )
        F = Assembler.assembleProjection!( i_traction_force_fc, interior_fc.fs_N_fc, interior_fc.index_fc, F )
    end

    weight_K( x, e, xi ) = chi( x ) * E( x ) * A( x ) * ( 1.0 / [ interior_fc.geom_field_fc.field_parametric_gradient( e, xi ) ] )
    i_K_fc = Integral.function_collection_integral( interior_fc.quadrature_fc, interior_fc.geom_field_fc, weight_K )
    K = Assembler.assembleInnerProduct!( i_K_fc, interior_fc.fs_dNds_fc, interior_fc.fs_dNds_fc, formulation.dis_strain_mat, interior_fc.index_fc, interior_fc.index_fc, K )

    weight_body_force( x, e, xi ) = interior_fc.f( x ) * chi( x ) * [ interior_fc.geom_field_fc.field_parametric_gradient( e, xi ) ]
    i_body_force_fc = Integral.function_collection_integral( interior_fc.quadrature_fc, interior_fc.geom_field_fc, weight_body_force )
    F = Assembler.assembleProjection!( i_body_force_fc, interior_fc.fs_N_fc, interior_fc.index_fc, F )

    return K, M, B, F, G, H
end

end
