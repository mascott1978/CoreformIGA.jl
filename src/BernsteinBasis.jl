module BernsteinBasis

using Plots

"""
    basisValue(i, p, xi)

hi
"""
function basisValue( i, p, xi )
    return binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
end

function basisValues( p, xi )
    return [ basisValue( i, p, xi ) for i in 0 : p ]
end

function basisDerivative( i, p, xi )
    return p * (basisValue( i - 1, p - 1, xi) - basisValue( i, p - 1, xi ) )
end

function basisDerivatives( p, xi )
    return [ basisDerivative( i, p, xi ) for i in 0 : p ]
end


function graph!( plt, i, p, f; subd = 10 )
    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = [ f( i, p, xi ) for xi in x ]
    Plots.plot!( plt, x, y )
    return plt
end

function graph!( plt, p, f; subd = 10 )
    x = [ xi for xi in LinRange( 0, 1, subd ) ];
    y = zeros( subd, p + 1 )
    count = 1
    for xi in x
        y[ count, : ] = f( p, xi );
        count += 1
    end
    for a = 1 : p + 1
        Plots.plot!( plt, x, y )
    end
    return plt
end

end

