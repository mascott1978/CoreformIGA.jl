module CoreformIGA

using Plots
using GaussQuadrature
using LinearAlgebra

include("BasisBernstein.jl")
include("BasisSpline.jl")
include("QuadratureGauss.jl")
include("Field.jl")
include("FlexRepresentationMethod1d.jl")
include("Examples.jl")

end
