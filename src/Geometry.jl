module Geometry

import LinearAlgebra
import ..NonlinearSolver

struct FunctionCollectionMapInversion
    geometric_map_inversion_predictor::Function
    norm::Function
    solver_linear::Function
    solver_nonlinear::Function
    geometric_map_inversion::Function
end

function function_collection_map_inversion_0d( bm_fc, field_fc )
    return FunctionCollectionMapInversion( x -> 1, 0.0,
                                           x -> 0.0,
                                           ( K, R ) -> 0.0,
                                           ( input, predictor, residual, tangent, solver_linear, norm, update ) -> 0.0,
                                           ( x, e_x, xi_x ) -> 1, 0.0 )
end

function function_collection_map_inversion_1d( bm_fc, field_fc )
    geometric_map_inversion_predictor = geometric_map_inversion_predictor_sequential_1d( bm_fc, field_fc )
    norm = LinearAlgebra.norm
    solver_linear( K, R ) = ( 1.0 / K ) * R
    solver_nonlinear = NonlinearSolver.newtonRaphsonIteration
    return FunctionCollectionMapInversion( geometric_map_inversion_predictor,
                                           norm,
                                           solver_linear,
                                           solver_nonlinear,
                                           geometric_map_inversion( bm_fc, field_fc, geometric_map_inversion_predictor,
                                                                    norm, solver_linear, solver_nonlinear ) )
end

function geometric_map_inversion_predictor_sequential_1d( bm_fc, field_fc )
    function geometric_map_inversion_predictor_sequential_1d( x )
        curr = bm_fc.element_count()
        for e in 1 : elem_n
            x_right = field_fc.field_value( e, 1.0 )
            if x < x_right
                curr = e
                break
            end
        end
        return curr, 0.5
    end
    return geometric_map_inversion_predictor_sequential_1d
end

function geometric_map_inversion( bm_fc,
                                  field_fc,
                                  geometric_map_inversion_predictor::Function,
                                  norm::Function,
                                  solver_linear::Function
                                  solver_nonlinear::Function )
    function geometric_map_inversion( x, e_x::Integer, xi_x )
        e::Integer, xi = geometric_map_inversion_predictor( x )
        residual( x, xi ) = x - field_fc.field_value( e, xi )
        tangent( xi ) = field_fc.field_parametric_gradient( e, xi )
        update( xi, delta_xi ) = xi + delta_xi
        return e, solver_nonlinear( x, x -> x, residual, tangent, solver_linear, norm, update )
    end
    return geometric_map_inversion
end

end
