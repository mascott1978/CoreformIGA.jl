module Evals

include("./basismesh.jl")

export buildEvals, evaluationpoint

#N and dN_ds are vectors of the spline evaluations at this point
struct evaluationpoint
       N
       dNds
       X
       dXds
       elem_id #element id number
       chi #indicator RV
       xi #gauss quadrature point
       w #weight
end


function storePoints( bm, quadratureFunc, eid, fict_penalty )
	 points = [ evaluationpoint( BasisMesh.computeN( bm, eid, quad_point.pt ),
	 	    		     BasisMesh.computedNds( bm, eid, quad_point.pt ),
	 	    		     BasisMesh.computeX( bm, eid, quad_point.pt ),
	 	    		     BasisMesh.computedXds( bm, eid, quad_point.pt ),
				     eid,
				     quad_point.chi ? 1 : fict_penalty,
                                     quad_point.pt,
                                     quad_point.wt )
     				     for quad_point in quadratureFunc( ( bm.pVec[eid]+1 ), eid ) ] #Fixme hardcoded p+1 quadrature
         return points
end

#builds evaluation points so we can assemply structures
function buildEvals( bm, quadratureFunc, fict_penalty )
	 evals = collect(Iterators.flatten([storePoints( bm, quadratureFunc, i, fict_penalty ) for i = 1:bm.elemN])) #get evaluation points for each element.
    return evals
end

end
