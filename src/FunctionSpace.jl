module FunctionSpace

struct FunctionCollectionFunctionSpace
    global_function_count_on_element::Function
    global_function_id_on_element::Function
    global_function_ids_on_element::Function
    global_basis_evaluator::Function
    global_basis_gradient_evaluator::Function
    geometric_map_inversion::Function
end

function function_collection_function_space( bm_fc, mi_fc, bs_fc )
    return FunctionCollectionFunctionSpace(  bm_fc.global_function_count_on_element,
                                             bm_fc.global_function_id_on_element,
                                             bm_fc.global_function_ids_on_element,
                                             bs_fc.global_basis_value,
                                             bs_fc.global_basis_parametric_gradient,
                                             mi_fc.geometric_map_inversion )
end

end
