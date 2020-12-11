using Plots

function computeBernstein(i,n,xi) return binomial(n,i)*xi^i*(1-xi)^(n-i)
end

function computeBersteinDeriv(i,n,xi)
    return binomial(n,i)*xi^(i-1)*(1-xi)^(-i+n-1)*(i-n*xi)
end

function getBernsteinVec(n, xi)
	return[computeBernstein(i,n,xi) for i = 0:n]
end

function getBernsteinDerivVec(n, xi)
   return[computeBersteinDeriv(i,n,xi) for i = 0:n]
end

function getSplineVec(p, xi, extraction_op)
	N = getBernsteinVec(p, xi)
	return extraction_op * N
end

function getSplineDerivVec(p, xi, extraction_op)
	dNds = getBernsteinDerivVec(p, xi)
	return extraction_op * dNds
end

function plotSplines( p, extraction_op, plt, range; steps=100 )
	 gr()
	 x = [i for i in LinRange(0,1,steps)]
	 xp = [i for i in LinRange(range[1],range[2],steps)]
	 y = [getSplineVec(p, xi, extraction_op) for xi in x]
	 yd = [getSplineDerivVec(p, xi, extraction_op) for xi in x]
	 #print( y )
	 for i=0:p
	     yp = [y[j][i+1] for j=1:100]
	     ydp = [yd[j][i+1] for j=1:100]
             a = i + 1
	     plt=plot!(xp,yp,label="B$a")
	     #plt=plot!(xp,ydp)
	 end
end

#Hardcoded for p2 c1
function getExtractionOp(p, elem_n, eid)
     if p != 2
     throw( "not implemented" )
     end
	 if elem_n == 1
	    return [ 1 0 0;
	    	     0 1 0;
		     0 0 1 ]
	 end
	 if eid == elem_n
	    return [ 1/2 0 0;
	    	     1/2 1 0;
		     0   0 1 ]
	 elseif eid == 1
	      return [ 1 0 0;
	      	       0 1 1/2;
		       0 0 1/2 ]
	 else
		return [ 1/2 0 0;
		         1/2 1 1/2;
			 0   0 1/2 ]
	 end
end
