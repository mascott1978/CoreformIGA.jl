module Integral

@enum Destination LHS_SYM=1 LHS_NONSYM=2 RHS=3

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
