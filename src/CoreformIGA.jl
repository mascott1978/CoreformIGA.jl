module CoreformIGA

include("1dCutCellDomains.jl")
include("Quadrature.jl")
include("BasisBernstein.jl")
include("BasisMesh.jl")
include("BEXTImport.jl")
include("BasisSpline.jl")
include("FunctionSpace.jl")
include("Field.jl")
include("Integral.jl")
include("NonlinearSolver.jl")
include("Assembler.jl")
include("Geometry.jl")
include("FlexRepresentationMethod.jl")
include("Viz.jl")
include("BoundaryCondition.jl")
include("./../element/Formulation1DSolid.jl")

end
