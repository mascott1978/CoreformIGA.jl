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
    return quadrature_point( i ) = layout.quadrature_points[ i ], layout.quadrature_weights[ i ]
end

function layout_gauss_legendre_0d()
    return Layout( [ 1.0 ], [ 1.0 ], [ 1 ] )
end

function layout_gauss_legendre_1d( element_count::Function, gauss_legendre_rules )
    quadrature_points = zeros( 0 )
    quadrature_weights= zeros( 0 )
    element_id_at_quadrature_point = zeros( 0 )
    for e in 1 : element_count()
        quadrature_points_e, quadrature_weights_e = gauss_legendre_rule[ e ]
        append!( quadrature_points, quadrature_points_e )
        append!( quadrature_weights, quadrature_weights_e )
        append!( element_id_at_quadrature_point, [ e for i in 1: length( quadrature_points_e ) ] )
    end
    return Layout( quadrature_points, quadratures_weights, element_id_at_quadrature_point )
end

function rule_gauss_legendre_1d( rule::Integer )
    qp, qw = GaussQuadrature.legendre( rule )
    return ( qp .+ 1.0 ) ./ 2.0, ( 1.0 / 2.0 ) .* qw
end

end

