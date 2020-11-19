module Assemble

function assembleInnerProduct!( geom_quadrature_rule,
                                geom_vals::Function,
                                geom_measure::Function,
                                weight::Function,
                                closest_parametric_point_test::Function,
                                basis_evaluator_test::Function,
                                basis_id_mapping_test::Function,
                                closest_parametric_point_trial::Function,
                                basis_evaluator_trial::Function,
                                basis_id_mapping_trial::Function, M )

    for i in geom_quadrature_rule
        x_i = geom_vals( i.e, i.qp )
        scaling_i = geom_measure( i.e, i.qp )

        e_test_i, xi_test_i = closest_parametric_point_test( x_i )
        M_i = basis_evaluator_test( e_test_i, xi_test_i )
        id_map_test_i = basis_id_mapping_test( e_test_i )

        e_trial_i, xi_trial_i = closest_parametric_point_trial( x_i )
        N_i = basis_evaluator_trial( e_trial_i, xi_trial_i )
        id_map_trial_i = basis_id_mapping_trial( e_trial_i )

        w = weight( x_i )
        for a in 1:size( M_i )[ 1 ]
            for b in 1:size( N_i )[ 1 ]
                M[ id_map_test_i[ a ], id_map_trial_i[ b ] ] += w * M_i[ a ] * N_i[ b ] * scaling_i * i.qw
            end
        end
    end
    return M

end

function assembleProjection!( geom_quadrature_rule,
                              geom_vals::Function,
                              geom_measure::Function,
                              weight::Function,
                              closest_parametric_point_test::Function,
                              basis_evaluator_test::Function,
                              basis_id_mapping_test::Function, F )

    for i in geom_quadrature_rule
        x_i = geom_vals( i.e, i.qp )
        scaling_i = geom_measure( i.e, i.qp )

        e_test_i, xi_test_i = closest_parametric_point_test( x_i )
        M_i = basis_evaluator_test( e_test_i, xi_test_i )
        id_map_test_i = basis_id_mapping_test( e_test_i )

        w = weight( x_i )
        for a in 1:size( M_i )[ 1 ]
            F[ id_map_test_i[ a ] ] += w * M_i[ a ] * scaling_i * i.qw
        end
    end
    return F
end

end

module NonlinearSolver

function newtonRaphsonIteration( input,
                                 predictor::Function,
                                 residual::Function,
                                 tangent::Function,
                                 solve::Function,
                                 norm::Function,
                                 update::Function;
                                 max_iter = 10,
                                 tol = 1e-6 )

    x = predictor( input )
    R = residual( input, x )
    iter = 0
    while norm( R ) > tol
        K = tangent( x )
        delta = solve( K, R )
        x = update( x, delta )
        R = residual( input, x )
        iter += 1
        if iter > max_iter
            throw( "Maximum iterations reached and no solution found!")
        end
    end
    return x
end

function uzawaIteration( K, M, B, F, G, H;
                         max_iter = 10 )

    d_curr = zeros( size( K )[ 1 ] + 1 )
    lambda_curr = zeros( 1 )
    for i in 1 : max_iter
        K = K + p_u .* M
        rhs_1 = ( F + p_u .* H - lambda_curr .* B  )
        d_curr = K \ rhs_1
        rhs_2 = p_u .* ( transpose(B) * d_curr - G)
        lambda_curr = lambda_curr + rhs_2
    end
    return d_curr, lambda_curr
end


end
