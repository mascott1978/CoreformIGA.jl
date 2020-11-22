module Viz

import Plots

function plotBasisLocal!( plt,
                          p::Int64,
                          local_basis::Function;
                          subd = 20,
                          domain = [ 0, 1 ] )

    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = zeros( subd, p + 1 )
    count = 1
    for xi in x
        y[ count, : ] = local_basis( p, xi );
        count += 1
    end
    for a in 1 : p + 1
        Plots.plot!( plt, ( domain[ 2 ] - domain[ 1 ] ) .* x .+ domain[ 1 ], y[ :, a ] )
    end
    return plt
end

function plotBasisSpline!( plt,
                           e::Int64,
                           p::Int64,
                           local_basis::Function,
                           extraction_operator::Function,
                           spline_basis::Function;
                           subd = 20,
                           domain = [ 0, 1 ] )

    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = zeros( subd, p + 1 )
    count = 1
    for xi in x
        y[ count, : ] = spline_basis( e, p, xi, extraction_operator, local_basis );
        count += 1
    end
    for a in 1 : p + 1
        Plots.plot!( plt, ( domain[ 2 ] - domain[ 1 ] ) .* x .+ domain[ 1 ], y[ :, a ] )
    end
    return plt
end

function plotField!( plt,
                     e::Int64,
                     p::Int64,
                     nodes_x::Function,
                     basis_evaluator_x::Function,
                     nodes_f::Function,
                     basis_evaluator_f::Function;
                     subd = 20,
                     domain = [ 0, 1 ] )

    xi = [ xi for xi in LinRange( 0, 1, subd ) ];
    x = [ evaluator( e, p, xi, nodes_x, basis_evaluator_x ) for xi in xi ]
    y = [ evaluator( e, p, xi, nodes_f, basis_evaluator_f ) for xi in xi ]
    Plots.plot!( plt, x, y )
    return plt
end

end
