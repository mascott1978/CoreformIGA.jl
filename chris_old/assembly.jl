module Assembly

export assembleStiffness, assembleConstraints, assembleForceVector, solveSystem, extractSolution

include("./basismesh.jl")
include("./inverse.jl")

import LinearAlgebra

function assembleStiffness( bm, evals )
	K = zeros( bm.funcN, bm.funcN )
	E = 1
	A = 1

	#loop through pre-evaluated quadrature points
	for pt in evals
	   dNdX = pt.dNds * (1/pt.dXds)
	   for ( ai, dNdX_a ) in enumerate( dNdX )
	       for( bi, dNdX_b ) in enumerate( dNdX )
	       	   K[ BasisMesh.getEqId( bm, pt.elem_id, ai ), BasisMesh.getEqId( bm, pt.elem_id, bi ) ] += pt.chi * E * A * dNdX_a * dNdX_b * pt.dXds * pt.w
	       end
	   end
	end
	return K
end

#constraints is a vector of pairs with [ x, value ]
function assembleConstraints( bm, constraints )
	n = size(constraints)[1]
	B = zeros( n, bm.funcN ) 
	FB = zeros( n )
	for (i, pt) in enumerate( constraints )
	    inverse_elem_id, s = Inverse.invertX( pt[1], bm )
	    N = BasisMesh.computeN( bm, inverse_elem_id, s )
	    for (b, N_b) in enumerate( N )
	         B[ i, BasisMesh.getEqId( bm, inverse_elem_id, b ) ] = N_b
	    end
	    FB[ i ] = pt[2]
	end
	return B, FB
end

#FIXME
function assembleForceVector( bm, evals, body_force, point_loads )
	 F = zeros( bm.funcN )
	 #assemble body force
	 #FIXME Make body force a function of X
	 for pt in evals
	     for ( ai, N_a ) in enumerate( pt.N )
	     	 F[ BasisMesh.getEqId( bm, pt.elem_id, ai ) ] += pt.chi * body_force * N_a * pt.dXds * pt.w 
	     end
	 end

	 #assemble traction
	 for pt in point_loads
	     inverse_elem_id, s = Inverse.invertX( pt[1], bm )
	     N = BasisMesh.computeN( bm, inverse_elem_id, s )
	     for ( ai, N_a ) in enumerate( N )
	     	 F[ BasisMesh.getEqId( bm, inverse_elem_id, ai ) ] += N_a * pt[2] 
	     end
	  end

	 return F
end

function solveSystem( K, B, FB, F )
    cN = size( B )[1]
    KB = [ [K; B] [B'; zeros( cN, cN )] ]
    Fs = [ F; FB ]
    print("KB ", KB, "\n")
    print("FS ", Fs, "\n")
    #print("KB ", size(KB), "\n")
    #print("FS ", size(Fs), "\n")
    d = KB\Fs
    return d, LinearAlgebra.cond( KB )
end

function extractSolution( d, bm )
	return d[1:bm.funcN]
end

end
