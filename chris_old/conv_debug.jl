#!/bin/bash
#=
exec julia -i --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#


include("./utils/quadrature.jl")
include("./utils/bernstein.jl")
include("./utils/basismesh.jl")
include("./utils/evals.jl")
include("./utils/assembly.jl")
include("./utils/inverse.jl")
include("./utils/ImportBEXT.jl")

constraints = [ [ 0.0, 0.0 ] ]#fix left end
#bodyload(x) = 0 #no body force (for now, but we will have one)
bodyload = 0 #no body force (for now, but we will have one)
point_loads = [ [ 1.0, 1.0e-3 ] ] #point load in x
physical_geometry = [0, 1]
fict_penalty = 1e-6

analytic_solution = 1e-3

filename = "/home/chris/cf/master/cf-iga-course/refinement_meshes/full_flex_bext_p2h1.json"
io = open(filename, "r")
file = read(io, String)

#import BasisMesh
bm = BasisMesh.buildMeshFromBEXT( file )
#print(bm)
evals = Evals.buildEvals( bm, Quadrature.returnFunctionExactGauss( bm, physical_geometry ), fict_penalty )
K = Assembly.assembleStiffness( bm, evals )
B, FB = Assembly.assembleConstraints( bm, constraints )
F = Assembly.assembleForceVector( bm, evals, bodyload, point_loads )
d, cond = Assembly.solveSystem( K, B, FB, F )
sol = Assembly.extractSolution( d, bm )
inverse_eid, inverse_s = Inverse.invertX( 1.0, bm ) # Find where the solution should be evaluated
tip_displacement = BasisMesh.computeField( bm, inverse_eid, inverse_s, sol )
error = abs( tip_displacement - analytic_solution )

println( "cond: ", cond )
print( "sol: ", sol, "\n")
print( error, "\n")


print("Press ENTER to exit")
readline()
#exit()
