module BasisSpline

struct FunctionCollection
    global_basis_value::Function
    global_basis_parametric_gradient::Function
end

function function_collection( bm_fc )
    return FunctionCollection( global_basis_value( bm_fc ),
                               global_basis_parametric_gradient( bm_fc ) )
end

function global_basis_value( bm_fc )
    return global_basis_value( e::Integer, xi::Array{ Float64, 1 } ) = global_basis_evaluator( e, xi, bm_fc.element_degree, bm_fc.local_function_count_on_element, bm_fc.extraction_operator_on_element, bm_fc.local_basis_value )
end

function global_basis_parametric_gradient( bm_fc )
    return global_basis_parametric_gradient( e::Integer, xi::Array{ Float64, 1 } ) = global_basis_evaluator( e, xi, bm_fc.element_degree, bm_fc.local_function_count_on_element, bm_fc.extraction_operator_on_element, bm_fc.local_basis_parametric_gradient )
end

function global_basis_evaluator( e::Integer, xi::Array{ Float64, 1 }, element_degree::Function, local_function_count_on_element::Function, extraction_operator_on_element::Function, local_basis_evaluator::Function )
    p = element_degree( e )
    B_xi = local_basis_evaluator( p, xi )
    C_e = extraction_operator_on_element( e )
    return C_e * B_xi
end

end
