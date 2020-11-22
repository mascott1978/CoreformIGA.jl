module FunctionSpace

struct FunctionCollectionFunctionSpace
    global_function_count_on_element::Function
    global_function_id_on_element::Function
    global_basis_evaluator::Function
    geometric_map_inversion::Function
end

function function_collection_function_space( bm_fc, mi_fc, global_basis_evaluator::Function )
    return FunctionCollectionFunctionSpace(  bm_fc.global_function_count_on_element,
                                             bm_fc.global_function_id_on_element,
                                             global_basis_evaluator,
                                             mi_fc.geometric_map_inversion )
end

end
