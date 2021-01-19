module Integral

struct FunctionCollectionIntegral
    quadrature_point_count::Function
    quadrature_point::Function
    integrand::Function
    destination::String
end

function function_collection_integral( q_fc, integrand, dest )
    return FunctionCollectionIntegral( q_fc.quadrature_point_count,
                                       q_fc.quadrature_point,
                                       integrand, dest )
end

end
