module Field

struct FunctionCollection
    node_on_element::Function
    nodes_on_element::Function
    nodes_on_domain::Function
    field_value::Function
    field_parametric_gradient::Function
end

function function_collection( bm_fc, bs_fc, nodes )

    return FunctionCollection( node_on_element( bm_fc, nodes ),
                               nodes_on_element( bm_fc, nodes ),
                               nodes_on_domain( nodes ),
                               field_value( bm_fc, bs_fc, nodes ),
                               field_parametric_gradient( bm_fc, bs_fc, nodes ) )
end

function node_on_element( bm_fc, nodes )
    return node_on_element( e, a ) = nodes[ bm_fc.global_function_id_on_element( e, a ) ]
end

function nodes_on_element( bm_fc, nodes )
    return nodes_on_element( e ) = [ nodes[ bm_fc.global_function_id_on_element( e, a ) ] for a in 1 : bm_fc.global_function_count_on_element( e ) ]
end

function nodes_on_domain( nodes )
    return nodes_on_domain( i ) = nodes[ :, i ]
end

function field_value( bm_fc, bs_fc, nodes )
    return field_value( e, xi ) = field_evaluator( e, xi, nodes_on_element( bm_fc, nodes ), bs_fc.global_basis_value )
end

function field_parametric_gradient( bm_fc, bs_fc, nodes )
    return field_parametric_gradient( e, xi ) = field_evaluator( e, xi, nodes_on_element( bm_fc, nodes ), bs_fc.global_basis_parametric_gradient )
end

function field_evaluator( e::Integer, xi::Array{ Float64, 1 }, nodes_on_element::Function, global_basis_evaluator::Function )
    coeffs = nodes_on_element( e )
    basis = global_basis_evaluator( e, xi )
    return sum( coeffs[ i ] * basis[ i ] for i in 1 : length( basis ) )
end

# function solveForNodes( node_n::Int64,
#                         param_quadrature_rule,
#                         param_vals::Function,
#                         param_measure::Function,
#                         input::Function,
#                         closest_parametric_point::Function,
#                         basis_evaluator::Function,
#                         basis_id_mapping::Function,
#                         assembler_lhs::Function,
#                         assembler_rhs::Function,
#                         solver::Function )
#     M = zeros( node_n, node_n )
#     F = zeros( node_n, 1 )
#     M = assembler_lhs( param_quadrature_rule,
#                        param_vals,
#                        param_measure,
#                        x -> 1.0,
#                        closest_parametric_point,
#                        basis_evaluator,
#                        basis_id_mapping,
#                        closest_parametric_point,
#                        basis_evaluator,
#                        basis_id_mapping, M )

#     F = assembler_rhs( param_quadrature_rule,
#                        param_vals,
#                        param_measure,
#                        input,
#                        closest_parametric_point,
#                        basis_evaluator,
#                        basis_id_mapping, F )

#     return solver( M, F )
# end


end
