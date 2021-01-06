module BasisBernstein

"""
    B(p, xi)

Return the value of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.
"""
function B( p::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    vec1 = B1d( p[ end ], xi[ end ] )

    if size( p, 1 ) == 1
        return vec1
    else
        vec2 = B( p[ 1:end-1 ], xi[ 1:end-1 ] )
        return reshape( kron( vec1, vec2 ), :, 1 )
    end
end

function B1d( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return reshape([ val( i, p, xi ) for i in 0 : p ],:,1)
end

"""
    dBdxi(p, xi)

Return the first derivative of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`B(p, xi)`](@ref)
"""
function dBdxi( p::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    r = zeros( 0 )

    function eval( pe::Int64, xie::Float64, ie, je )
        if ie == je
            return dBdxi1d( pe, xie )
        else
            return B1d( pe, xie )
        end
    end

    for i in 0 : length( p ) - 1
        vec2 = [ 1 ]
        for j in 0 : length( p ) - 1
            vec1 = eval( p[ end - j ], xi[ end - j ], i, j )
            vec2 = kron( vec1, vec2 )
        end
        prepend!( r, vec2 )
    end
    return reshape( r, :, length( p ) )
end

function dBdxi1d( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return reshape([ p * (val( i - 1, p - 1, xi) - val( i, p - 1, xi ) ) for i in 0 : p ],:,1)
end


"""
    dBdxi(p, n, xi)

Return the nth derivative of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`dBdxi(p, xi)`](@ref)
"""
function dBdxi( p::Array{ Int64, 1 }, n::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    vec2 = [ 1 ]
    for i in 0 : length( p ) - 1
        vec1 = dBdxi1d( p[ end - i ], n[ end - i ], xi[ end - i ] )
        vec2 = kron( vec1, vec2 )
    end
    return reshape( vec2, :, length( p ) )
end

function dBdxi1d( p::Int64, n::Int64, xi::Float64 )
    function val( i, p, n, xi )
        if n == 0
            if i > p || i < 0
                return 0;
            else
                return binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
            end
        end

        if p < 0
            return 0
        end

        return p*( val( i - 1, p - 1, n - 1, xi ) - val( i, p - 1, n - 1, xi ) )
    end

    return reshape([ val( i, p, n, xi ) for i in 0 : p ],:,1)
end
end
