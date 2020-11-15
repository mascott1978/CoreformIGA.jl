module CoreformIGA

using Plots
using GaussQuadrature
using LinearAlgebra

include("BasisBernstein.jl")
include("BasisBspline_UniformHMaxK.jl")
include("QuadratureGauss.jl")
include("Field.jl")
include("FlexRepresentationMethod1d.jl")

function test1()
    deg = 2
    elem_n = 20
    quad_rule = 3
    quad_rules = [ quad_rule for i in 1:elem_n ]
    cad_domain = [ 1, 2 ]
    flex_domain = [ 0, 3 ]
    p_cad = 1e-12
    p_u = 1
    E = 1
    A = 1
    function load( x )
        return 0
    end
    function traction( x )
        return 1
    end
    function constraint( x )
        return 0
    end
    FlexRepresentationMethod1d.solve( deg, elem_n, quad_rules, flex_domain,
                                      p_cad, cad_domain, E, A, load, traction, p_u, constraint )
end


end
