
module BasisSpline

using Plots
using ..BasisBernstein

struct Layout
    domain
    lengths
    degrees
    smoothnesses
    ops
    EG
    func_n
end

function buildUniformHMaxK( p, elem_n; domain = [0, 1] )
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
    return Layout( domain, lengths, degrees, smoothnesses, ops, EG, func_n )
end

function N( layout, e, xi )
    p = layout.degrees[ e ]
    B_xi = BasisBernstein.B( p, xi )
    C_e = layout.ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ]
    return C_e * B_xi
end

function dNdxi( layout, e, xi )
    p = size( layout.ops )[ 1 ] - 1
    dBdxi_xi = BasisBernstein.dBdxi( p, xi )
    C_e = layout.ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ]
    return C_e * dBdxi_xi
end

function plot!( plt, layout, e, f; subd = 20, domain = [ 0, 1 ] )
    p = size( layout.ops )[ 1 ] - 1
    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = zeros( subd, p + 1 )
    count = 1
    for xi in x
        y[ count, : ] = f( layout, e, xi );
        count += 1
    end
    for a = 1 : p + 1
        Plots.plot!( plt, ( domain[ 2 ] - domain[ 1 ] ) .* x .+ domain[ 1 ], y[ :, a ] )
    end
    return plt
end

function nodesEquallySpaced( layout )
    return [ node for node in layout.domain[ 1 ]:(layout.domain[ 2 ] - layout.domain[ 1 ]) / ( size(layout.lengths )[ 1 ] + 1):layout.domain[ 2 ] ]
end

end

