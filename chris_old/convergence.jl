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


h_max = 6
k_max = 6

output_data = [[] for i=1:k_max]
output_cond = [[] for i=1:k_max]
output_h = [ 2^i for i=1:6 ]

#Linear mappings r=constant
for k=1:6 #k refinement
    output_k = zeros( h_max )
    k_cond = zeros( h_max )
    for h=1:h_max #h refinement
        filename = "/home/chris/cf/master/cf-iga-course/refinement_meshes/full_flex_bext_p$(k)h$(h).json"
	print(filename, "\n")
        #print(filename,"\n")
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
        output_k[h] = abs( tip_displacement - analytic_solution )
        k_cond[h] = cond
    end
    output_data[k] = output_k
    output_cond[k] = k_cond
end

print(output_h)
print("\n")
for i=1:k_max
    print("         p = ", i, " ", output_data[i], "\n")
    print("cond for p = ", i, " ", output_cond[i], "\n")
end

plt=plot(xaxis=:log, yaxis=:log)
for i=1:k_max
    plt=plot!(output_h, output_data[i], label=string(i) )
end
gui(plt)


#constraints = [ [ 0.0, 0 ], [1, 0.5], [3, 0.5], [4, 1] ]
#fict_penalty = 1e-9
#physical_geometry = [0, 4]
#flex_geometry = [ 0, 4 ]
#basis_mesh = BasisMesh.buildUniformMeshGreville( 2, 4, flex_geometry )
#
#evals = Evals.buildEvals( basis_mesh, Quadrature.returnFunctionGauss( basis_mesh, physical_geometry ), fict_penalty )
#
#K = Assembly.assembleStiffness( basis_mesh, evals )
#B, FB = Assembly.assembleConstraints( basis_mesh, constraints )
#F = Assembly.assembleForceVector( basis_mesh, 0, [] )
#
#d = Assembly.solveSystem( K, B, FB, F )
#
##print(d)
#s = Assembly.extractSolution( d, basis_mesh )
#plt=plot()
##print( "\ns", s, "\n")
#plt=plot( title="Solution" )
#BasisMesh.plotField( basis_mesh, plt, s, name="Displacement" )
#gui(plt)

print("Press ENTER to exit")
readline()
exit()
