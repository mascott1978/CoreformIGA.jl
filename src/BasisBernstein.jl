module BasisBernstein

"""
    B(p, xi)

Return the value of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.
"""
function B( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return [ val( i, p, xi ) for i in 0 : p ]
end

"""
    dBdxi(p, xi)

Return the first derivative of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`B(p, xi)`](@ref)
"""
function dBdxi( p::Int64, xi::Float64 )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return [ p * (val( i - 1, p - 1, xi) - val( i, p - 1, xi ) ) for i in 0 : p ]
end

end
