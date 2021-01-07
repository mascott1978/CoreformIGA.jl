module BasisBernstein

function tensor_product( components )
    vec1 = components[ end ]()

    if length( components ) == 1
        return vec1
    else
        vec2 = tensor_product( components[ 1:end-1 ] )
        return reshape( kron( vec1, vec2 ), :, 1 )
    end
end

"""
    B(p, xi)

Return the value of the tensor product Bernstein basis of degree ``p`` evaluated at ``\\xi``.
This can be multidimensional and so both ``p`` and ``n`` are arrays to specify the basis order in each direction and the number of derivatives in each direction, respectively.
"""
function B( p::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    components = [ ()->B1d( p[ i ], xi[ i ] ) for i in 1:length( p ) ]
    return tensor_product( components )
end

"""
    B1d(p, xi)

Return the value of the 1d Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.
"""
function B1d( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return reshape([ val( i, p, xi ) for i in 0 : p ],:,1)
end



"""
    dBdxi(p,xi)

Return the gradient of the tensor product Bernstein basis of degree ``p`` evaluated at ``\\xi``.
This can be multidimensional and so both ``p`` and ``n`` are arrays to specify the basis order in each direction and the number of derivatives in each direction, respectively.
"""
function dBdxi( p::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    r = zeros( 0 )
    for i in 1 : length( p )
        components = [ ( i == j ) ? ( ()->dBdxi1d( p[ i ], xi[ i ] ) ) : ( ()->B1d( p[ j ], xi[ j ] ) ) for j in 1:length( p ) ]
        append!( r, tensor_product( components ) )
    end
    return reshape( r, :, length( p ) )
end

"""
    dBdxi1d(p, xi)

Return the first derivative of the 1d Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`B1d(p, xi)`](@ref)
"""
function dBdxi1d( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return reshape([ p * (val( i - 1, p - 1, xi) - val( i, p - 1, xi ) ) for i in 0 : p ],:,1)
end



"""
    dnBdxin(p, n, xi)

Return the ``nth`` derivative of the tensor product Bernstein basis of degree ``p`` evaluated at ``\\xi``.
This can be multidimensional and so both ``p`` and ``n`` are arrays to specify the basis order in each direction and the number of derivatives in each direction, respectively.
"""
function dnBdxin( p::Array{ Int64, 1 }, n::Array{ Int64, 1 }, xi::Array{ Float64, 1 } )
    components = [ ()->dnBdxin1d( p[ i ], n[ i ], xi[ i ] ) for i in 1:length( p ) ]
    return tensor_product( components )
end

"""
    dnBdxin1d(p, n, xi)

Return the nth derivative of the 1d Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`dBdxi1d(p, xi)`](@ref)
"""
function dnBdxin1d( p::Int64, n::Int64, xi::Float64 )
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
