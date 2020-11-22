module CoreformIGA

include("NonlinearSolver.jl")
include("Assembler.jl")

include("Quadrature.jl")

include("BasisBernstein.jl")
include("BasisMesh.jl")
include("BasisSpline.jl")

include("FunctionSpace.jl")
include("Field.jl")
include("Geometry.jl")
include("Integral.jl")

include("FlexRepresentationMethod1d.jl")
include("Viz.jl")

end
