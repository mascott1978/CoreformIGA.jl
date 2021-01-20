module Quadrature

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

# 0d code
function layout_gauss_legendre_0d()
    return Layout( [ [ 1.0 ] ], [ 1.0 ], [ 1 ] )
end

# 1d code
function layout_gauss_legendre_1d( element_count::Function, element_degree::Function, inverse_map, d_bc_fc, t_bc_fc, #=What's the point of this second input=#gauss_legendre_rules )
    quadrature_points = Any[]
    quadrature_weights= zeros( 0 )
    element_id_at_quadrature_point = zeros( Int64, 0 )

    domain_process_eval( e, x0, x1 ) = begin
        quadrature_points_e, quadrature_weights_e = rule_gauss_legendre_1d( element_degree( e )[ 1 ] + 1, x0, x1 )
        append!( quadrature_points, quadrature_points_e )
        append!( quadrature_weights, quadrature_weights_e )
        append!( element_id_at_quadrature_point, [ e for i in 1: length( quadrature_points_e ) ] )
    end
    cuts = layout_cut_cell_domains1d( element_count, d_bc_fc, t_bc_fc, inverse_map )
    iterate_domains( cuts, domain_process_eval )
    return Layout( quadrature_points, quadrature_weights, element_id_at_quadrature_point )
end

function rule_gauss_legendre_1d( rule::Integer, x0::Float64=0.0, x1::Float64=1.0 )
    qp, qw = GaussQuadrature.legendre( rule )
    return [ [ x0 * ( 1.0 - qp[ i ] ) / 2.0 + x1 * ( qp[ i ] + 1.0 ) / 2.0 ] for i in 1 : length( qp ) ], ( ( x1 - x0 ) / 2.0 ) .* qw
end

struct ElementCuts
    cells::Array{ Array{ Float64, 1 }, 1 }
end

function iterate_domains( cuts::ElementCuts, process_domain::Function )
    for e in 1 : size( cuts.cells, 1 )
        for cut in 1 : ( size( cuts.cells[ e ], 1 ) - 1 )
            x_0 = cuts.cells[ e ][ cut ]
            x_1 = cuts.cells[ e ][ cut + 1 ]
            process_domain( e, x_0, x_1 )
        end
    end
end

function layout_cut_cell_domains1d( element_count::Function, d_bc_fc, n_bc_fc, inverse_map )
    cells::Array{ Array{ Float64, 1 }, 1 } = fill([], element_count())
    for e in 1:size( cells, 1 )
        cells[ e ] = [ 0.0, 1.0 ]
    end

    nodes_d = d_bc_fc
    # for n in 1 : size( nodes_d, 1 )
    if ~isempty( nodes_d )
        e, x = inverse_map.geometric_map_inversion( [ nodes_d, 0, 0 ], 0, 0.0 )
        if ( abs( x[ 1 ] - 0.0 ) > 1e-10 ) && ( abs( x[ 1 ] - 1.0 ) > 1e-10 )
            append!( cells[ e ], x[ 1 ] )
        end
    end

    nodes_n =  n_bc_fc 
    # for n in 1 : size( nodes_n, 1 )
    if ~isempty( nodes_n )
        e, x = inverse_map.geometric_map_inversion( [ nodes_n, 0, 0 ], 0, 0.0 )
        if ( abs( x[ 1 ] - 0.0 ) > 1e-10 ) && ( abs( x[ 1 ] - 1.0 ) > 1e-10 )
            append!( cells[ e ], x[ 1 ] )
        end
    end

    for e in 1 : size( cells, 1 )
        sort!( cells[ e ] )
    end
    return ElementCuts( cells )
end

end

