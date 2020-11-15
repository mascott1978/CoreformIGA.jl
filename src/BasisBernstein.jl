module BasisBernstein

using Plots


"""
    B(p, xi)

Return the value of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.
"""
function B( p, xi )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return [ val( i, p, xi ) for i in 0 : p ]
end

"""
    dBdxi(p, xi)

Return the first derivative of the Bernstein basis of degree ``p`` evaluated at ``\\xi \\in [0,1]``.

See also: [`B(p, xi)`](@ref)
"""
function dBdxi( p, xi )
    val( i, p, xi ) = binomial( p, i ) * xi^i * ( 1 - xi )^( p - i )
    return [ p * (val( i - 1, p - 1, xi) - val( i, p - 1, xi ) ) for i in 0 : p ]
end

"""
    plot!(plt, p, f; subd=10)

Plot the values or first derivatives (as indicated by the passed in function [`B(p, xi)`](@ref) or [`dBdxi(p, xi)`](@ref)) of the Bernstein basis of degree ``p`` over ``[0,1]``.

See also: [`B(p, xi)`](@ref), [`dBdxi(p, xi)`](@ref)
"""
function plot!( plt, p, f; subd = 10 )
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
