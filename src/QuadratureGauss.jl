module QuadratureGauss

struct quad_point
    e::Int
    qp::Float64
    qw::Float64
end

using GaussQuadrature

function legendreRule( rule )
    qp, qw = GaussQuadrature.legendre( rule )
    return ( qp .+ 1.0 ) ./ 2.0, ( 1.0 / 2.0 ) .* qw
end

function legendreLayout( rules )
    qp = zeros( 0 )
    qw= zeros( 0 )
    IE = zeros( 0 )
    e = 1
    for rule in rules
        qp_e, qw_e = legendreRule( rule )
        append!( qp, qp_e )
        append!( qw, qw_e )
        append!( IE, [e for i in 1:rule] )
        e += 1
    end
    quad_points = Array{ quad_point, 1 }( undef, size( qp )[ 1 ] )
    for i in 1:size( qp )[ 1 ]
        quad_points[ i ] = quad_point( IE[ i ], qp[ i ], qw[ i ] )
    end
    return quad_points
end

end

