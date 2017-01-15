"""Module for integration of polynomials over 3D volumes and surfaces"""
using LAR


function M(alpha, beta)
    a = 0
    for l=1:(alpha + 2)
        a += binomial(alpha+1,l) * (-1)^l/(l+beta+1)
    end
    return a/(alpha + 1)
end

""" The main integration routine """
function TT(tau::Array{Float64,2}, alpha, beta, gamma, signedInt=false)
	vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
	a = va - vo
	b = vb - vo
	s1 = 0.0
	for h=0:alpha
		for k=0:beta
			for m=0:gamma
				s2 = 0.0
				for i=0:h 
					s3 = 0.0
					for j=0:k
						s4 = 0.0
						for l=0:m
							s4 += binomial(m,l) * a[3]^(m-l) * b[3]^l * M( 
								h+k+m-i-j-l, i+j+l )
						end
						s3 += binomial(k,j) * a[2]^(k-j) * b[2]^j * s4
					end
					s2 += binomial(h,i) * a[1]^(h-i) * b[1]^i * s3;
				end
				s1 += binomial(alpha,h) * binomial(beta,k) * binomial(gamma,m) * 			
						vo[1]^(alpha-h) * vo[2]^(beta-k) * vo[3]^(gamma-m) * s2
			end
		end
	end
	c = cross(a,b)
	if signedInt == true
		return s1 * vecnorm(c) * sign(c[3])
	else
		return s1 * vecnorm(c)
	end	
end

""" Basic integration functions """
function II(P, alpha, beta, gamma, signedInt=false)
    w = 0
    V, FV = P
    if typeof(P) == PyCall.PyObject
    	if typeof(V) == Array{Any,2}
    		V = V'
    	end
    	if typeof(FV) == Array{Any,2}
    		FV = [FV[k,:] for k=1:size(FV,1)]
    		FV = FV+1
    	end
    end
    if typeof(FV) == Array{Int64,2}
    	FV = [FV[:,k] for k=1:size(FV,2)]
    end
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        if size(tau,2) == 3
        	term = TT(tau, alpha, beta, gamma, signedInt)
        	if signedInt
        		w += term
        	else
        		w += abs(term)
        	end
        elseif size(tau,2) > 3
        	println("ERROR: FV[$(i)] is not a triangle")
        else
        	println("ERROR: FV[$(i)] is degenerate")
        end
    end    
    return w
end

function III(P, alpha, beta, gamma)
    w = 0
    V, FV = P
    if typeof(P) == PyCall.PyObject
    	if typeof(V) == Array{Any,2}
    		V = V'
    	end
    	if typeof(FV) == Array{Any,2}
    		FV = [FV[k,:] for k=1:size(FV,1)]
    		FV = FV+1
    	end
    end
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
        a = va - vo
        b = vb - vo
        c = cross(a,b)
        w += c[1]/vecnorm(c) * TT(tau, alpha+1, beta, gamma)
    end
    return w/(alpha + 1)
end


""" Surface and volume integrals """
function Surface(P, signedInt=false)
    return II(P, 0, 0, 0, signedInt)
end

function Volume(P)
    return III(P, 0, 0, 0)
end

#v,(vv,ev,fv,cv) = p.larCuboids((1,1,1),true)
#V = hcat([Array{Float64,1}(v[k,:]) for k=1:size(v,1)]...)
#FV = hcat([Array{Int64,1}(fv[k,:]+1) for k=1:size(fv,1)]...)
#EV = hcat([Array{Int64,1}(ev[k,:]+1) for k=1:size(ev,1)]...)
#model1 = Any[V,FV,EV]
#P = V,[FV[:,k] for k=1:size(FV,2)]
#surface(P,false)
#
#
#""" Surface integration """
#function surfIntegration(larModel)
#    V,FV,EV = model
#    FE = crossRelation(FV,EV)
#    if typeof(FV) == Array{Int64,2}
#    	FV = [FV[:,k] for k=1:size(FV,2)]
#    end
#    if typeof(V) == Array{Int64,2}
#    	if size(V,1) == 2
#    		V = vcat(V,zeros(1,size(V,2)))
#    	end
#    end
#    cochain = []
#    triangles = []
#    faceVertPairs = []
#	for face=1:length(FE)
#		push!(faceVertPairs, hcat([EV[:,e] for e in FE[face]]...))
#		row = [faceVertPairs[face][1] for k=1:length(FE[face])]
#		push!(triangles, vcat(faceVertPairs[face],row'))
#        P = V,triangles[face]
#        area = surface(P,false) 
#        push!(cochain,area)
#    end
#    return cochain
#end
    
#    TODO: after having implemented ∂_3/∂_2
#def signedIntSurfIntegration(model,signedInt=False)
#    V,FV,EV = model
#    V = [v+[0.0] if len(v)==2 else v for v in V]
#    cochain = []
#    for face in FV
#        triangles = AA(C(AL)(face[0]))(TRANS([face[1-1],face[2]]))
#        P = V,triangles
#        area = Surface(P,signedInt) 
#        cochain += [area]
#    return cochain
#


#""" Terms of the Euler tensor """
#def FirstMoment(P)
#    out = [None]*3
#    out[0] = III(P, 1, 0, 0)
#    out[1] = III(P, 0, 1, 0)
#    out[2] = III(P, 0, 0, 1)
#    return out
#
#def SecondMoment(P)
#    out = [None]*3
#    out[0] = III(P, 2, 0, 0)
#    out[1] = III(P, 0, 2, 0)
#    out[2] = III(P, 0, 0, 2)
#    return out
#
#def InertiaProduct(P)
#    out = [None]*3
#    out[0] = III(P, 0, 1, 1)
#    out[1] = III(P, 1, 0, 1)
#    out[2] = III(P, 1, 1, 0)
#    return out
#
#""" Vectors and covectors of mechanical interest """
#def Centroid(P)
#    out = [None]*3
#    firstMoment = FirstMoment(P)
#    volume = Volume(P)
#    out[0] = firstMoment[0]/volume
#    out[1] = firstMoment[1]/volume
#    out[2] = firstMoment[2]/volume
#    return out
#
#def InertiaMoment(P)
#    out = [None]*3
#    secondMoment = SecondMoment(P)
#    out[0] = secondMoment[1] + secondMoment[2]
#    out[1] = secondMoment[2] + secondMoment[0]
#    out[2] = secondMoment[0] + secondMoment[1]
#    return out

