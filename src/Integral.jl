module Integral

struct FunctionCollectionIntegral
    quadrature_point_count::Function
    quadrature_point::Function
    geometry_value::Function
    integral_weight::Function
end

function function_collection_integral( q_fc, geom_fc, integral_weight )
    return FunctionCollectionIntegral( q_fc.quadrature_point_count,
                                       q_fc.quadrature_point,
                                       geom_fc.field_value,
                                       integral_weight )
end

end
