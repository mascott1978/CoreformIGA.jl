#!/bin/bash
#=
exec julia -i --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#
#@show ARGS


include("./utils/quadrature.jl")
include("./utils/bernstein.jl")
include("./utils/basismesh.jl")
include("./utils/evals.jl")
#include("./utils/geom.jl")
include("./utils/assembly.jl")
include("./utils/inverse.jl")
include("./utils/ImportBEXT.jl")

#constraints at certain points along the bar [position, constraint]
#constraints = [ [ 0.0, 0 ], [1, 0.5], [3, 0.5], [4, 1] ]
constraints = [ [ 0.0, 0 ], [1, 0.5] ]
#used for indicator R.V.
fict_penalty = 1e-9
#domain of physical/CAD geometry
physical_geometry = [0, 1]
#domain of flex geometry (has to contain physical geometry)
flex_geometry = [ 0, 4 ]

#builds the mesh using graville points unifrom continuity and degree k=p-1
basis_mesh = BasisMesh.buildUniformMeshGreville( 2, 7, flex_geometry )

#returnFunctionGauss returns a function that will return an array containing the quadrature point, weight, and
#indicator of whether the quadrature point is in the domain (physical)
evals = Evals.buildEvals( basis_mesh, Quadrature.returnFunctionExactGauss( basis_mesh, physical_geometry ), fict_penalty )

K = Assembly.assembleStiffness( basis_mesh, evals )
B, FB = Assembly.assembleConstraints( basis_mesh, constraints )
F = Assembly.assembleForceVector( basis_mesh, evals, 0, [] )

d = Assembly.solveSystem( K, B, FB, F )

#print(d)
s = Assembly.extractSolution( d, basis_mesh )

inverse_eid, inverse_s = Inverse.invertX( 1.0, basis_mesh ) # Find where the solution should be evaluated
tip_displacement = BasisMesh.computeField( basis_mesh, inverse_eid, inverse_s, s )
print(tip_displacement)
plt=plot()
#print( "\ns", s, "\n")
plt=plot( title="Solution" )
BasisMesh.plotField( basis_mesh, plt, s, name="Displacement" )
gui(plt)

print("Press ENTER to exit")
readline()
exit()
