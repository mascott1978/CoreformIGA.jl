module Field

function localize( e::Int64,
                   e_func_n::Int64,
                   nodes::Function )

    return [ nodes( e, a ) for a in 1 : e_func_n ]
end

function evaluator( e::Int64,
                    e_func_n::Int64,
                    xi::Float64,
                    nodes::Function,
                    basis_evaluator::Function,
                    inner_product::Function )

    N_e = basis_evaluator( e, xi )
    P_e = localize( e, e_func_n, nodes )
    return inner_product( P_e, N_e )
end

end

module Geometry

function nodesFromParameterization( node_n::Int64,
                                    param_quadrature_rule,
                                    param_vals::Function,
                                    param_measure::Function,
                                    parameterization::Function,
                                    closest_parametric_point::Function,
                                    basis_evaluator::Function,
                                    basis_id_mapping::Function,
                                    assembler_lhs::Function,
                                    assembler_rhs::Function,
                                    solver::Function )
    M = zeros( node_n, node_n )
    F = zeros( node_n, 1 )
    M = assembler_lhs( param_quadrature_rule,
                       param_vals,
                       param_measure,
                       x -> 1.0,
                       closest_parametric_point,
                       basis_evaluator,
                       basis_id_mapping,
                       closest_parametric_point,
                       basis_evaluator,
                       basis_id_mapping, M )

    F = assembler_rhs( param_quadrature_rule,
                       param_vals,
                       param_measure,
                       parameterization,
                       closest_parametric_point,
                       basis_evaluator,
                       basis_id_mapping, F )

    return solver( M, F )
end

end
