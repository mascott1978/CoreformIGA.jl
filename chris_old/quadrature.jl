module Quadrature

include("./basismesh.jl")
include("./inverse.jl")
using GaussQuadrature

export QuadraturePoint, returnFunctionGauss, returnFunctionExactGauss, legendreRule

struct QuadraturePoint
    pt #quadrature point (coordinate)
    wt #quadrature weight
    chi #True means inside, False means outside -> indicator R.V.
end

function legendreRule( order )
    qp, qw = GaussQuadrature.legendre( order )
    return ( qp .+ 1.0 ) ./ 2.0, ( 1.0 / 2.0 ) .* qw
end

#Physical geom is a pair [ L_0, L_1 ]
#this function returns the function GaussChi. GaussChi returns an array of instances of the QuadraturePoint
#struct
function returnFunctionGauss( bm, physical_geom )
    function GaussChi( order, eid )
        xi, w = legendreRule( order )
        return [ QuadraturePoint( xi[i], w[i], Inverse.Chi( bm, physical_geom, eid, xi[i] )  ) for i=1:order]
    end
    return GaussChi
end

function returnFunctionExactGauss( bm, physical_geom )
    function GaussChi( order, eid )
        segments = Inverse.getSegments( bm, eid, physical_geom ) #returns array of pairs (xi_0, xi_1)
        pts = []
        xi, w = legendreRule( order )
        for segment in segments
            seg_pts, seg_wts = changeUnitInterval( xi, w, segment )
            append!( pts, [ QuadraturePoint( seg_pts[i], seg_wts[i], Inverse.Chi( bm, physical_geom, eid, seg_pts[i] ) )  for i=1:order ] )
        end
        return pts
    end
    return GaussChi
end


function changeUnitInterval( qpts, qwts, ab )
    b_minus_a = ab[2] - ab[1];
    #perform a change of variables, and add the scaling x' to the wts
    for pt in qpts
        pt = pt * b_minus_a + ab[1]
    end
    qwts = qwts * b_minus_a
    return qpts, qwts
end

end
