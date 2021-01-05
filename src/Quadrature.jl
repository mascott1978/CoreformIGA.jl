module Quadrature

using ..CutCellDomains1d

using GaussQuadrature

struct FunctionCollectionQuadrature
    quadrature_point_count::Function
    quadrature_point::Function
end

struct Layout
    quadrature_points
    quadrature_weights
    element_id_at_quadrature_point
end

function function_collection_quadrature( layout::Layout )

    return FunctionCollectionQuadrature( quadrature_point_count( layout ),
                                         quadrature_point( layout ) )
end

function quadrature_point_count( layout::Layout )
    return quadrature_point_count() = length( layout.quadrature_points )
end

function quadrature_point( layout::Layout )
    return quadrature_point( i ) = layout.element_id_at_quadrature_point[ i ], layout.quadrature_points[ i ], layout.quadrature_weights[ i ]
end

function layout_gauss_legendre_0d()
    return Layout( [ [ 1.0 ] ], [ 1.0 ], [ 1 ] )
end

function layout_gauss_legendre_1d( element_count::Function, element_degree::Function, inverse_map, d_bc_fc, t_bc_fc, #=What's the point of this second input=#gauss_legendre_rules )
    quadrature_points = Any[]
    quadrature_weights= zeros( 0 )
    element_id_at_quadrature_point = zeros( Int64, 0 )

    domain_process_eval( e, x0, x1 ) = domain_process( element_degree, quadrature_points, quadrature_weights, element_id_at_quadrature_point, e, x0, x1 )
    layout = CutCellDomains1d.layout_cut_cell_domains( element_count, d_bc_fc, t_bc_fc, inverse_map )
    CutCellDomains1d.iterate_domains( layout, domain_process_eval )
    return Layout( quadrature_points, quadrature_weights, element_id_at_quadrature_point )
end

# function rule_gauss_legendre_1d( rule::Integer )
#     qp, qw = GaussQuadrature.legendre( rule )
#     return ( qp .+ 1.0 ) ./ 2.0, ( 1.0 / 2.0 ) .* qw
# end

function rule_gauss_legendre_1d( rule::Integer, x0::Float64=0.0, x1::Float64=1.0 )
    qp, qw = GaussQuadrature.legendre( rule )
    return [ [ x0 * ( 1.0 - qp[ i ] ) / 2.0 + x1 * ( qp[ i ] + 1.0 ) / 2.0 ] for i in 1 : length( qp ) ], ( ( x1 - x0 ) / 2.0 ) .* qw
end

function domain_process( element_degree::Function, quadrature_points, quadrature_weights, element_id_at_quadrature_point, e, x0, x1 )
    quadrature_points_e, quadrature_weights_e = rule_gauss_legendre_1d( element_degree( e )[ 1 ] + 1, x0, x1 )
    append!( quadrature_points, quadrature_points_e )
    append!( quadrature_weights, quadrature_weights_e )
    append!( element_id_at_quadrature_point, [ e for i in 1: length( quadrature_points_e ) ] )
end

end

