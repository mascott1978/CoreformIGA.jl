module BasisMesh

struct Layout
    dim
    domain
    starts
    lengths
    degrees
    smoothnesses
    ops
    EG
    func_n
end

function buildBspline0d( val = 1.0,
                         domain = [ 0 ] )

    return Layout( 0, domain, [ domain[ 1 ] ], [ 0 ], [ 0 ], [ -1 ], [ val ], [ 1 ], 1 )
end

function buildBspline1d( p::Int64,
                         elem_n::Int64;
                         domain = [ 0, 1 ] )

    function op( e, p, elem_n )
        if p == 0
            return 1;
        elseif p == 1
            return [ 1 0;
                     0 1 ]
        elseif p == 2
            if elem_n == 1
                return [ 1 0 0;
                         0 1 0;
                         0 0 1 ]
            end
            if e == elem_n
                return [ 1/2 0 0;
                         1/2 1 0;
                         0   0 1 ]
            elseif e == 1
                return [ 1 0 0;
                         0 1 1/2;
                         0 0 1/2 ]
            else
                return [ 1/2 0 0;
                         1/2 1 1/2;
                         0   0 1/2 ]
            end
        elseif p == 3
            if elem_n == 1
                return [ 1 0 0 0;
                         0 1 0 0;
                         0 0 1 0;
                         0 0 0 1 ]
            elseif elem_n == 2
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 1/2;
                             0 0 0 1/4 ]
                else
                    return [ 1/4 0 0 0;
                             1/2 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            elseif elem_n == 3
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 7/12;
                             0 0 0 1/6 ]
                elseif e == 2
                    return [ 1/4 0 0 0;
                             7/12 2/3 1/3 1/6;
                             1/6 1/3 2/3 7/12;
                             0 0 0 1/4 ]
                else
                    return [ 1/6 0 0 0;
                             7/12 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            else
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 7/12;
                             0 0 0 1/6 ]
                elseif e == 2
                    return [ 1/4 0 0 0;
                             7/12 2/3 1/3 1/6;
                             1/6 1/3 2/3 2/3;
                             0 0 0 1/6 ]
                elseif e > 2 && e < elem_n - 1
                    return [ 1/6 0 0 0;
                             2/3 2/3 1/3 1/6;
                             1/6 1/3 2/3 2/3;
                             0 0 0 1/6 ]
                elseif e == elem_n - 1
                    return [ 1/6 0 0 0;
                             2/3 2/3 1/3 1/6;
                             1/6 1/3 2/3 7/12;
                             0 0 0 1/4 ]
                else
                    return [ 1/6 0 0 0;
                             7/12 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            end
        else
            error( "invalid degree specified in extraction operator construction.")
        end
    end
    elem_size = ( domain[ 2 ] - domain[ 1 ] ) / elem_n
    lengths = [ elem_size for i in 1 : elem_n ]
    starts = zeros( elem_n );
    start = domain[ 1 ]
    for e in 1 : elem_n
        starts[ e ] = start
        start += lengths[ e ]
    end

    degrees = [ p for i in 1 : elem_n ]
    smoothnesses = [ p - 1 for i in 1 : elem_n ]
    smoothnesses[ 1 ] = -1
    append!( smoothnesses, -1 )
    ops = zeros( p + 1, ( p + 1 ) * elem_n );
    EG = zeros( Int32, elem_n, p + 1 )
    for e in 1 : elem_n
        ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ] = op( e, p, elem_n )
        EG[ e, : ] = e : e + p
    end
    func_n = elem_n + p
    return Layout( 1, domain, starts, lengths, degrees, smoothnesses, ops, EG, func_n )
end
end


module BasisSpline

function evaluate( e::Int64,
                   p::Int64,
                   xi::Float64,
                   extraction_operator::Function,
                   local_basis_vals::Function )

    B_xi = local_basis_vals( p, xi )
    C_e = extraction_operator( e )
    return C_e * B_xi
end

end

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
