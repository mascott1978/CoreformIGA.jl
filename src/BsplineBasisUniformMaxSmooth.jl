module BsplineBasisUniformMaxSmooth

include("BernsteinBasis.jl")
using Plots

function build( p, elem_n  )
    ops = zeros( p + 1, ( p + 1 ) * elem_n );
    EG = zeros( Int32, elem_n, p + 1 )
    for e = 1 : elem_n
        ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ] = op( e, p, elem_n )
        EG[ e, : ] = e : e + p
    end
    return ops, EG
end

function basisValue( e, a, xi, ops )
    return basisValues( e, xi, ops )[ a, : ]
end

function basisValues( e, xi, ops )
    p = size( ops )[ 1 ] - 1
    B = BernsteinBasis.basisValues( p, xi )
    C_e = ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ]
    return C_e * B
end

function basisDerivative( e, a, xi, ops )
    return basisDerivatives( e, xi, ops )[ a, : ]
end

function basisDerivatives( e, xi, ops )
    p = size( ops )[ 1 ] - 1
    B = BernsteinBasis.basisDerivatives( p, xi )
    C_e = ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ]
    return C_e * B
end

function graph!( plt, e, a, f, ops; subd = 10, flex_domain = [ 0, 1 ] )
    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = [ f( e, a, xi, ops ) for xi in x ]
    Plots.plot!( plt, ( flex_domain[ 2 ] - flex_domain[ 1 ] ) .* x .+ flex_domain[ 1 ], y )
    return plt
end

function graph!( plt, e, f, ops; subd = 20, flex_domain = [ 0, 1 ] )
    p = size( ops )[ 1 ] - 1
    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = zeros( subd, p + 1 )
    count = 1
    for xi in x
        y[ count, : ] = f( e, xi, ops );
        count += 1
    end
    for a = 1 : p + 1
        Plots.plot!( plt, ( flex_domain[ 2 ] - flex_domain[ 1 ] ) .* x .+ flex_domain[ 1 ], y[ :, a ] )
    end
    return plt
end

function knotVector( p, elem_n, flex_domain )
    kv = [ knot for knot in flex_domain[ 1 ]:(flex_domain[2]-flex_domain[1])/elem_n:flex_domain[2] ]
    prepend!( kv, [ flex_domain[ 1 ] for i in 1 : p - 1 ] )
    append!( kv, [ flex_domain[ 2 ] for i in 1 : p - 1 ] )
    return kv
end

function grevilleNodes( p, elem_n, flex_domain )
    kv = knotVector( p, elem_n, flex_domain )
    nodes = [ sum( kv[ i : i + p - 1 ]) / p for i = 1 : elem_n + p ]
    return nodes
end

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

end

