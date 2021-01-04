module NonlinearSolver

function newtonRaphsonIteration( input,
                                 predictor::Function,
                                 residual::Function,
                                 tangent::Function,
                                 solver_linear::Function,
                                 norm::Function,
                                 update::Function;
                                 max_iter = 10,
                                 tol = 1e-6 )

    curr = predictor( input ) # second filter
    R = residual( input, curr )
    iter = 0
    while norm( R ) > tol
        K = tangent( curr )
        delta = solver_linear( K, R )
        curr = update( curr, delta )
        R = residual( input, curr )
        iter += 1
        if iter > max_iter
            throw( "Maximum iterations reached and no solution found!")
        end
    end
    return curr
end

function uzawaIteration( K, M, B, F, G, H;
                         max_iter = 10, p_u = 1 )

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
