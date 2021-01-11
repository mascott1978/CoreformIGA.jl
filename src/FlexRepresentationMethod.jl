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
                   dirichlet_bcs,
                   neumann_bcs,
                   formulation,
                   chi::Function,
                   penalty_constraint::Function,
                   E::Function,
                   A::Function )

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
