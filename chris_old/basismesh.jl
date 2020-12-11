module BasisMesh

include("./bernstein.jl")
include("./ImportBEXT.jl")

export basismesh, buildMeshFromBEXT, buildUniformMesh, getNodes, getEqId, computeN, computedNds, computeX, computedXds, plotBasis, plotField
using Plots

mutable struct basismesh
    elemN #The number of elements in the mesh
    pVec #The degree of each element. Is of length elemN
    kVec #The continuity between elements. Is length elemN-1
    extraction #An elemN vec that stores the dense extraction operator for each element
    funcN #Total number of functions in this mesh
    nodes #The given geometric mapping
    indexing #an elemN vec that stores the global function index for each local funtion (in order)
end

#Read the basis and mesh from an externally generated spline
function buildMeshFromBEXT( json )
    p, ee, en, ns = ImportBEXT.bext( json )
    elemN = size(en)[1]
    #print( "p ", p, "\n" )
    pVec = [x[1] for x in p]
    kVec = [] #This is empty for BEXT
    extraction = [ vcat(transpose.(vecs)...) for vecs in ee ]
    #print("extraction", extraction, "\n")
    funcN = size(ns)[1]
    nodes = [x[1] for x in ns]
    indexing = [x .+ 1 for x in en]

    #return 0

    return basismesh( elemN, pVec, kVec, extraction, funcN, nodes, indexing )
end

#Construct the basis and mesh for a uniform p, k=p-1 mesh with a control points located at greville points
function buildUniformMeshGreville( p, elemN, flex_geom )
    num_funcs = elemN+p #global function number
    nodes = zeros( num_funcs ) #initial array with 0s for nodes (evaluation points)
    indexing = [ [ i+j for j=0:p ] for i=1:elemN ] #local to global function index map
    extraction = [getExtractionOp(p, elemN, i) for i=1:elemN] #get all extraction operators
    #builds greville points for nodes/control points/evaluatoin points
    #linrange is the same as np.linspace() in python(I think)
    if p == 1
	nodes = [x for x in LinRange(flex_geom[1], flex_geom[2], elemN + 1)]
    end
    dist = flex_geom[2] - flex_geom[1] #length of flex geometry

    if p == 2
	for i = 1:elemN+p
	    if i == 1
	       	nodes[i] = flex_geom[1]
	    elseif i == 2
	       	nodes[i] = ( 0.5 / elemN ) * dist + flex_geom[1]
	    elseif i == elemN+p-1
	       	nodes[i] = ( ( elemN - 0.5 ) / elemN ) * dist + flex_geom[1]
	    elseif i == elemN+p
	       	nodes[i] = flex_geom[2]
	    else
	        nodes[i] = ( ( i-2 + 0.5 ) / elemN ) * dist + flex_geom[1]
	    end
	end
    else throw( "Not Implemented" )
    end
    #return basismesh instance
    return basismesh( elemN, [p for i=1:elemN], [p-1 for i=1:elemN-1], extraction, num_funcs, nodes, indexing )
end


#bm = basismesh
# BASIS MESH ROUTINES

#returns the nodes for a given element
function getNodes( bm, eid )
    index = bm.indexing[eid]
    return [bm.nodes[i] for i in index]
end

#gets global function id (for a local function index)
function getEqId( bm, elem_id, local_func_index )
    return bm.indexing[elem_id][local_func_index]
end

#s is a parametric location
#evaluates spline functions for an elment of degree p at location s (using extraction operator)
function computeN( bm, eid, s )
    p = bm.pVec[eid]
    return getSplineVec( p, s, bm.extraction[eid] )
end

#evaluaes spline derivatives for an elment of degree p at location s (using extraction operator)
function computedNds( bm, eid, s )
    p = bm.pVec[eid]
    return getSplineDerivVec( p, s, bm.extraction[eid] )
end

function computeX( bm, eid, s )
    N = computeN( bm, eid, s )
    nodes = getNodes( bm, eid )
    #print( "eid ", eid, "\n")
    #print( "N ", N, "\n")
    #print( "nodes ", nodes, "\n")
    return nodes' * N
end

function computedXds( bm, eid, s )
    dNds = computedNds( bm, eid, s )
    nodes = getNodes( bm, eid )
    return nodes' * dNds
end

function computeField( bm, eid, s, field )
    N = computeN( bm, eid, s )
    index = bm.indexing[eid]
    field_nodes = [field[i] for i in index]
    field_nodes = vcat( transpose.(field_nodes)... )
    return field_nodes' * N
end

function plotBasis( bm, plt; steps=100 )
    gr()
    vals = [[] for i=1:bm.funcN*2]
    for i=1:bm.elemN
	r_min = i - 1
	r_max = i
	x = [j for j in LinRange(0,1,steps)]
	xp = [j for j in LinRange(r_min,r_max,steps)]
	y = [computeN( bm, i, xi) for xi in x]
	for j=1:bm.pVec[i]+1
	    yp = [y[k][j] for k=1:steps]
	    fid = getEqId( bm, i, j )
	    append!( vals[ fid * 2 ], yp )
	    append!( vals[ fid * 2 - 1 ], xp )
        end
    end
    #print(vals)
    for i=1:bm.funcN
	plt=plot!(vals[i*2-1],vals[i*2], label="N$i")
    end
end

function plotGeometry2d( bm, plt; steps=100 )
    gr()
    x = []
    field_val = []
    for i=1:bm.elemN
	append!( field_val, [computeField( bm, i, xi, bm.nodes ) for xi in LinRange(0,1,steps)] )
    end
    plt=plot!([i[1] for i in field_val], [i[2] for i in field_val], label="Geometry")
end

function plotField( bm, plt, field; steps=100, name="field" )
    gr()
    x = []
    field_val = []
    for i=1:bm.elemN
	append!( x, [computeX( bm, i, xi ) for xi in LinRange(0,1,steps)] )
	append!( field_val, [computeField( bm, i, xi, field ) for xi in LinRange(0,1,steps)] )
    end
    #print(x, "\n")
    #print(field_val, "\n")
    plt=plot!(x, field_val, label=name)
end

end
