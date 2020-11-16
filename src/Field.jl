module Field

using LinearAlgebra
using Plots
using ..BasisSpline

function localize( layout, nodes, e )
    p = layout.degrees[ e ]
    f_e = [ nodes[ layout.EG[ e, a ] ] for a in 1 : p + 1 ]
    return f_e
end

function F( layout, nodes, e, xi )
    N_e = BasisSpline.N( layout, e, xi )
    P_e = localize( layout, nodes, e )
    return LinearAlgebra.dot( P_e, N_e )
end

function dFdxi( layout, nodes, e, xi )
    dNdxi_e = BasisSpline.dNdxi( layout, e, xi )
    P_e = localize( layout, nodes, e )
    return LinearAlgebra.dot( P_e, dNdxi_e )
end

function plot!( plt, layout, nodes_x, nodes_f, e, f; subd = 10 )
    xi = [ xi for xi in LinRange( 0, 1, subd ) ];
    x = [ F( layout, nodes_x, e, xi ) for xi in xi ]
    y = [ f( layout, nodes_f, e, xi ) for xi in xi ]
    Plots.plot!( plt, x, y )
    return plt
end

function closestPoint( layout, nodes, x; max_iter = 10, tol = 1e-6 )
    elem_n = size( layout.EG )[ 1 ]
    inverse_element_id = elem_n
    for e = 1 : elem_n
        x_s1 = F( layout, nodes, e, 1 )
        if x < x_s1
            inverse_element_id = e
            break
        end
    end

    s_i = 0.5
    R = x - F( layout, nodes, inverse_element_id, s_i )
    iter = 0
    while abs(R) > tol
        dXds = dFdxi( layout, nodes, inverse_element_id, s_i )
        delta_s = R / dXds
        s_i += delta_s
        R = x - F( layout, nodes, inverse_element_id, s_i )
        iter += 1
        if iter > max_iter
            throw( "Maximum iterations reached and no solution found")
        end
    end
    return inverse_element_id, s_i
end

end

