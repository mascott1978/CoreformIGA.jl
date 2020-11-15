module Field

using LinearAlgebra
using Plots
using ..BasisBspline_UniformHMaxK

function localizeFieldToElement( e, nodes, EG )
    p = size( EG )[ 2 ] - 1
    f_e = [ nodes[ EG[ e, a ] ] for a in 1 : p + 1 ]
end

function dXdxi( e, xi, nodes, ops, EG )
    dNdxi_e = BasisBspline_UniformHMaxK.basisDerivatives( e, xi, ops )
    P_e = localizeFieldToElement( e, nodes, EG )
    return LinearAlgebra.dot( P_e, dNdxi_e )
end

function X( e, xi, nodes, ops, EG )
    N_e = BasisBspline_UniformHMaxK.basisValues( e, xi, ops )
    P_e = localizeFieldToElement( e, nodes, EG )
    return LinearAlgebra.dot( P_e, N_e )
end

function graph!( plt, e, nodes_x, nodes_y, ops, EG; subd = 10 )
    xi = [ xi for xi in LinRange( 0, 1, subd ) ];
    x = [ X( e, xi, nodes_x, ops, EG ) for xi in xi ]
    y = [ X( e, xi, nodes_y, ops, EG ) for xi in xi ]
    Plots.plot!( plt, x, y )
    return plt
end

function closestPoint( x, nodes, ops, EG; max_iter = 10, tol = 1e-6 )
    elem_n = size( EG )[ 1 ]
    inverse_element_id = elem_n
    for e = 1 : elem_n
        x_s1 = X( e, 1, nodes, ops, EG )
        if x < x_s1
            inverse_element_id = e
            break
        end
    end

    s_i = 0.5
    R = x - X( inverse_element_id, s_i, nodes, ops, EG )
    iter = 0
    while abs(R) > tol
        dXds = dXdxi( inverse_element_id, s_i, nodes, ops, EG )
        delta_s = R / dXds
        s_i += delta_s
        R = x - X( inverse_element_id, s_i, nodes, ops, EG )
        iter += 1
        if iter > max_iter
            throw( "Maximum iterations reached and no solution found")
        end
    end
    return inverse_element_id, s_i
end

end

