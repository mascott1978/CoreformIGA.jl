module Integral

"""
# Specify the destination of the integration in ``KU = F``
- `LHS_SYM`: the integration results in a symmetric matrix which needs to be assembled into ``K`` once
- `LHS_NONSYM`: the integration results in a non-symmetric matrix which needs to be assembled into ``K`` twice
- `RHS`: the integration results in a column vector which needs to be assembled into ``F``.
"""
@enum Destination LHS_SYM=1 LHS_NONSYM=2 RHS=3

"""
# Members:
- `quadrature_point_count`: get the number of quadrature points
- `quadrature_point`: get the quadrature points and weights
- `integrand`: get the integrand.
- `destination`: get the destination of the integration in ``KU = F``
"""
struct FunctionCollectionIntegral
    quadrature_point_count::Function
    quadrature_point::Function
    integrand::Function
    destination::Destination
end

function function_collection_integral( q_fc, integrand, dest )
    return FunctionCollectionIntegral( q_fc.quadrature_point_count,
                                       q_fc.quadrature_point,
                                       integrand, dest )
end

end
