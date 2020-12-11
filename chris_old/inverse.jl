module Inverse

export invertX, getSegments, Chi

include("./basismesh.jl")

#FIXME assumes that elements are ordered from left to right
function invertX( x, bm; tol=1e-9, max_iter=20, compare_tol = 1e-9 )

    print("Inverting ", x, "\n")
    #filter
    #find the correct element id by checking s = 1 for each element
    inverse_element_id = bm.elemN + 1 #INDEX OUT OF RANGE
    #print( "starting element: ", inverse_element_id, "\n")
    for i=1:bm.elemN
	x_s1 = BasisMesh.computeX( bm, i, 1.0 ) + compare_tol
	x_s0 = BasisMesh.computeX( bm, i, 0.0 ) - compare_tol
        #print( x_s1, " ", x_s0, "\n" )
        #print( "x_s1 ", x_s1, "\n" )
        #print( i, "\n" )
	if x < x_s1 && x > x_s0
	    inverse_element_id = i
	    break
	end
        #print("On Element: ", inverse_element_id, "\n")
    end
    print("On Element: ", inverse_element_id, " of ", bm.elemN, "\n")

    #set initial guess
    s_i = 0.5
    R = x - BasisMesh.computeX( bm, inverse_element_id, s_i )
    iter = 0
    print( "Starting Residual ", R, "\n" )
    while abs(R) > tol
	dXds = BasisMesh.computedXds( bm, inverse_element_id, s_i )
	delta_s = R / dXds
	s_i += delta_s
	R = x - BasisMesh.computeX( bm, inverse_element_id, s_i )
        print( "Current Residual ", R, "\n" )
	iter += 1
	if iter > max_iter
	    print( "current s_i", s_i )
	    throw( "Maximum iterations reached and no solution found")
	end
    end
    return inverse_element_id, s_i
end

function Chi( bm, physical_geom, eid, pt )
    X = BasisMesh.computeX( bm, eid, pt )
    if X < physical_geom[1] || X > physical_geom[2]
        return false
    else
        return true
    end
end

function getCutX( bm, eid, physical_geom, tol = 1e-9 )
    left_end = BasisMesh.computeX( bm, eid, 0.0 )
    right_end = BasisMesh.computeX( bm, eid, 1.0 )
    #print( "\nleft_end ", left_end, "\n")
    #print( "\nright_end ", right_end, "\n")
    #print( "\ngeom", physical_geom, "\n")
    if left_end < physical_geom[1] + tol && right_end > physical_geom[1] - tol
        return physical_geom[1]
    elseif left_end < physical_geom[2] + tol && right_end > physical_geom[2] - tol
        return physical_geom[2]
    else
        throw( "ERROR no cut element" )
    end
end

# NOTE Hardcoded xi from 0 to 1 NOTE
function getSegments( bm, eid, physical_geom )
    left_end_inside = Chi( bm, physical_geom, eid, 0.0 )
    right_end_inside = Chi( bm, physical_geom, eid, 1.0 )

    if left_end_inside == right_end_inside
        return [[0.0, 1.0]]
    else
        cut_elem, cut_point_xi = invertX( getCutX( bm, eid, physical_geom ), bm )
        return [[0.0, cut_point_xi], [cut_point_xi, 1.0]]
    end
end

end
