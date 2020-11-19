module FlexRepresentationMethod1d

using LinearAlgebra
using Plots
using ..BasisBernstein
using ..BasisSpline
using ..QuadratureGauss
using ..Field
using ..NonlinearSolver
using ..Assemble

function processBasis0d( layout )

    quadrature_rule = [ QuadratureGauss.quad_point( 1, 0, 1.0 ) ]
    extraction_operator( e ) = layout.ops
    spline_basis_N( e, xi ) = [ layout.ops[ 1 ] ]
    spline_basis_dNdxi( e, xi ) = 1
    id_map( e ) = [ 1 ]
    return quadrature_rule, extraction_operator, spline_basis_N, spline_basis_dNdxi, id_map
end

function processGeometry0d( nodes )

    node_map( e, a ) = nodes[ 1 ]
    geom_x( e, xi ) = nodes[ 1 ]
    geom_dxdxi( e, xi ) = 1
    closest_parametric_point( input ) = [ 1, 0 ]
    return node_map, geom_x, geom_dxdxi, closest_parametric_point
end

function processBasis1d( layout,
                         quad_rules )

    quadrature_rule = QuadratureGauss.legendreLayout( quad_rules )
    extraction_operator( e ) = layout.ops[ :, ( e - 1 ) * ( layout.degrees[ e ] + 1 ) + 1 : e * ( layout.degrees[ e ] + 1 ) ]
    spline_basis_N( e, xi ) = BasisSpline.evaluate( e, layout.degrees[ e ], xi,
                                                    extraction_operator,
                                                    BasisBernstein.B )
    spline_basis_dNdxi( e, xi ) = BasisSpline.evaluate( e, layout.degrees[ e ], xi,
                                                        extraction_operator,
                                                        BasisBernstein.dBdxi )
    id_map( e ) = [ layout.EG[ e, a ] for a in 1 : layout.degrees[ e ] + 1 ]
    return quadrature_rule, extraction_operator, spline_basis_N, spline_basis_dNdxi, id_map
end

function processGeometry1d( layout,
                            nodes,
                            spline_basis_N::Function,
                            spline_basis_dNdxi::Function )

    node_map( e, a ) = nodes[ layout.EG[ e, a ] ]
    geom_x( e, xi ) = Field.evaluator( e, layout.degrees[ e ] + 1, xi, node_map, spline_basis_N, LinearAlgebra.dot )
    geom_dxdxi( e, xi ) = Field.evaluator( e, layout.degrees[ e ] + 1, xi, node_map, spline_basis_dNdxi, LinearAlgebra.dot )
    function predictor( input )
        elem_n = size( layout.degrees )[ 1 ] 
        curr = elem_n
        for e in 1 : elem_n
            x_right = geom_x( e, 1.0 )
            if input < x_right
                curr = e
                break
            end
        end
        return Union{Int64, Float64}[ curr, 0.5 ]
    end
    residual( input, e_and_xi ) = input - geom_x( e_and_xi[ 1 ], e_and_xi[ 2 ] )
    tangent( e_and_xi ) = geom_dxdxi( e_and_xi[ 1 ], e_and_xi[ 2 ] )
    solve( K, R ) = ( 1.0 / K ) * R
    update( e_and_xi, delta_xi ) = Union{Int64, Float64}[ e_and_xi[ 1 ], e_and_xi[ 2 ] + delta_xi ]
    norm( R ) = abs( R )
    closest_parametric_point( input ) = NonlinearSolver.newtonRaphsonIteration( input, predictor, residual,
                                                                                tangent, solve, norm, update )
    return node_map, geom_x, geom_dxdxi, closest_parametric_point
end

function assemble( layout_interior,
                   nodes_interior,
                   quad_rules_interior,
                   layout_constraint_bdry,
                   nodes_constraint_bdry,
                   layout_traction_bdry,
                   nodes_traction_bdry,
                   chi::Function,
                   penalty_constraint::Function,
                   constraint::Function,
                   E::Function,
                   A::Function,
                   load::Function,
                   traction::Function )

    quadrature_rule_interior,
    extraction_operator_interior,
    spline_basis_N_interior,
    spline_basis_dNdxi_interior,
    id_map_interior = processBasis1d( layout_interior,
                                      quad_rules_interior )

    node_map_interior,
    geom_x_interior,
    geom_dxdxi_interior,
    closest_parametric_point_interior = processGeometry1d( layout_interior, nodes_interior,
                                                           spline_basis_N_interior,
                                                           spline_basis_dNdxi_interior )

    function spline_basis_dNdx( e, xi )
        dNdxi = spline_basis_dNdxi_interior( e, xi )
        dxdxi = geom_dxdxi_interior( e, xi )
        return dNdxi * ( 1.0 / dxdxi )
    end

    quadrature_rule_constraint_bdry,
    extraction_operator_constraint_bdry,
    spline_basis_N_constraint_bdry,
    spline_basis_dNdxi_constraint_bdry,
    id_map_constraint_bdry = processBasis0d( layout_constraint_bdry )

    node_map_constraint_bdry,
    geom_x_constraint_bdry,
    geom_dxdxi_constraint_bdry,
    closest_parametric_point_constraint_bdry = processGeometry0d( nodes_constraint_bdry )

    quadrature_rule_traction_bdry,
    extraction_operator_traction_bdry,
    spline_basis_N_traction_bdry,
    spline_basis_dNdx_traction_bdry,
    id_map_traction_bdry = processBasis0d( layout_traction_bdry )

    node_map_traction_bdry,
    geom_x_traction_bdry,
    geom_dxdxi_traction_bdry,
    closest_parametric_point_traction_bdry = processGeometry0d( nodes_traction_bdry )

    weight_K( x ) = chi( x ) * E( x ) * A( x )
    weight_M( x ) = penalty_constraint( x )
    weight_B( x ) = weight_M( x )
    weight_body_force( x ) = load( x ) * chi( x )
    weight_traction_force( x ) = traction( x )
    weight_G( x ) = constraint( x )
    weight_H( x ) = penalty_constraint( x ) * constraint( x )

    u_func_n = layout_interior.func_n
    lambda_func_n = layout_constraint_bdry.func_n

    K = zeros( u_func_n, u_func_n )
    M = zeros( u_func_n, u_func_n )
    B = zeros( lambda_func_n, u_func_n )
    F = zeros( u_func_n, 1 )
    G = zeros( lambda_func_n, 1 )
    H = zeros( u_func_n, 1 )

    K = Assemble.assembleInnerProduct!( quadrature_rule_interior,
                                        geom_x_interior,
                                        geom_dxdxi_interior,
                                        weight_K,
                                        closest_parametric_point_interior,
                                        spline_basis_dNdx_interior,
                                        id_map_interior,
                                        closest_parametric_point_interior,
                                        spline_basis_dNdx_interior,
                                        id_map_interior, K )

    M = Assemble.assembleInnerProduct!( quadrature_rule_constraint_bdry,
                                        geom_x_constraint_bdry,
                                        geom_dxdxi_constraint_bdry,
                                        weight_M,
                                        closest_parametric_point_interior,
                                        spline_basis_N_interior,
                                        id_map_interior,
                                        closest_parametric_point_interior,
                                        spline_basis_N_interior,
                                        id_map_interior, M )

    B = Assemble.assembleInnerProduct!( quadrature_rule_constraint_bdry,
                                        geom_x_constraint_bdry,
                                        geom_dxdxi_constraint_bdry,
                                        weight_B,
                                        closest_parametric_point_constraint_bdry,
                                        spline_basis_N_constraint_bdry,
                                        id_map_constraint_bdry,
                                        closest_parametric_point_interior,
                                        spline_basis_N_interior,
                                        id_map_interior, B )

    F = Assemble.assembleProjection!( quadrature_rule_interior,
                                      geom_x_interior,
                                      geom_dxdxi_interior,
                                      weight_body_force,
                                      closest_parametric_point_interior,
                                      spline_basis_N_interior,
                                      id_map_interior, F )

    F = Assemble.assembleProjection!( quadrature_rule_traction_bdry,
                                      geom_x_traction_bdry,
                                      geom_dxdxi_traction_bdry,
                                      weight_traction_force,
                                      closest_parametric_point_interior,
                                      spline_basis_N_interior,
                                      id_map_interior, F )

    G = Assemble.assembleProjection!( quadrature_rule_constraint_bdry,
                                      geom_x_constraint_bdry,
                                      geom_dxdxi_constraint_bdry,
                                      weight_G,
                                      closest_parametric_point_constraint_bdry,
                                      spline_basis_N_constraint_bdry,
                                      id_map_constraint_bdry, G )

    H = Assemble.assembleProjection!( quadrature_rule_constraint_bdry,
                                      geom_x_constraint_bdry,
                                      geom_dxdxi_constraint_bdry,
                                      weight_H,
                                      closest_parametric_point_interior,
                                      spline_basis_N_interior,
                                      id_map_interior, F )

    return K, M, B, F, G, H
end

end

# function assembleK!( quad_points, dNdxi, X, dXdxi, chi, id_map, E, A, K )
#     for i in quad_points
#         dNdxi_i = dNdxi( i.e, i.qp )
#         X_i = X( i.e, i.qp )
#         dXdxi_i = dXdxi( i.e, i.qp )
#         p_cad = chi( X_i )
#         dNdx_i = dNdxi_i .* ( 1.0 / dXdxi_i )
#         for a in 1:size( dNdx_i )[ 1 ]
#             for b in 1:size( dNdx_i )[ 1 ]
#                 res = p_cad * E * A * dNdx_i[ a ] * dNdx_i[ b ] * dXdxi_i * i.qw
#                 K[ id_map( i.e, a ), id_map( i.e, b ) ] += p_cad * E * A * dNdx_i[ a ] * dNdx_i[ b ] * dXdxi_i * i.qw
#             end
#         end
#     end
#     return K
# end

# function assembleM!( closest_point, cad_domain, N, id_map, p_u, M )
#     e_left, xi_left = closest_point( cad_domain[ 1 ] )
#     N_left = N( e_left, xi_left )
#     for a in 1:size( N_left )[ 1 ]
#         for b in 1:size( N_left )[ 1 ]
#             M[ id_map( e_left, a ), id_map( e_left, b ) ] += p_u * N_left[ a ] * N_left[ b ]
#         end
#     end
#     return M
# end

# function assembleB!( closest_point, cad_domain, N, id_map, B )
#     e_left, xi_left = closest_point( cad_domain[ 1 ] )
#     N_left = N( e_left, xi_left )
#     for a in 1:size( N_left )[ 1 ]
#         B[ id_map( e_left, a ) ] += N_left[ a ]
#     end
#     return B
# end

# function assembleF!( quad_points, closest_point, cad_domain, N, X, dXdxi, chi, id_map, load, traction, F )
#     for i in quad_points
#         N_i = N( i.e, i.qp )
#         X_i = X( i.e, i.qp )
#         dXdxi_i = dXdxi( i.e, i.qp )
#         p_cad = chi( X_i )
#         for a in 1:size( N_i )[ 1 ]
#             F[ id_map( i.e, a ) ] += p_cad * N_i[ a ] * load( X_i ) * dXdxi_i * i.qw
#         end
#     end
#     e_right, xi_right = closest_point( cad_domain[ 2 ] )
#     N_right = N( e_right, xi_right )
#     X_right = X( e_right, xi_right )
#     for a in 1:size( N_right )[ 1 ]
#         F[ id_map( e_right, a ) ] += N_right[ a ] * traction( X_right )
#     end
#     return F
# end

# function assembleG!( closest_point, cad_domain, X, constraint, G )
#     e_left, xi_left = closest_point( cad_domain[ 1 ] )
#     X_left = X( e_left, xi_left )
#     G[ 1 ] = constraint( X_left )
#     return G
# end

# function assembleH!( closest_point, cad_domain, N, X, id_map, p_u, constraint, H )
#     e_left, xi_left = closest_point( cad_domain[ 1 ] )
#     N_left = N( e_left, xi_left )
#     X_left = X( e_left, xi_left )
#     for a in 1:size( N_left )[ 1 ]
#         H[ id_map( e_left, a ) ] += p_u * N_left[ a ] * constraint( X_left )
#     end
#     return H
# end
