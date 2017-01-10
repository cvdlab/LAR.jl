#include("merge.jl")
#include("boundary.jl")

""" The `spacePartition` function takes as input a **non-valid** (with the meaning used in solid modeling field --- see~\cite{Requicha:1980:RRS:356827.356833}) `LAR` model of dimension d-1, i.e.~a triple `(V,FV,EV)`, The `spacePartition` function returns the **valid** `LAR` boundary model `(W,FW,EW)` of the space partition induced by `FV`. First an array `buckets` indexed on faces is computed. It contains the subsets of faces with greatest probability of intersecting each indexing face, respectively. 
"""

""" Half-line crossing test """
function crossingTest(new,old,count,status)
    if status == 0
        status = new
        count += 0.5
    else
        if status == old
        	count += 0.5
        else
        	count -= 0.5
        end
        status = 0
    end
end

function setTile(box)
	tiles = [[9,1,5],[8,0,4],[10,2,6]]
	b1,b2,b3,b4 = box
	function tileCode(point)
		x,y = point
		code = 0
		if y>b1 code=code|1 end
		if y<b2 code=code|2 end
		if x>b3 code=code|4 end
		if x<b4 code=code|8 end
		return code
	end
	return tileCode
end

# Point in polygon classification #
function pointInPolygonClassification(V,EV)
    function pointInPolygonClassification0(pnt)
        x,y = pnt
        xmin,xmax,ymin,ymax = x,x,y,y
        tilecode = setTile([ymax,ymin,xmax,xmin])
        count,status = 0,0
    
        for (k,edge) in enumerate(EV)
            p1,p2 = V[:,edge[1]],V[:,edge[2]]
            (x1,y1),(x2,y2) = p1,p2
            c1,c2 = tilecode(p1),tilecode(p2)
            c_edge, c_un, c_int = c1$c2, c1|c2, c1&c2
            
            if (c_edge == 0) & (c_un == 0) return "p_on" 
            elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
            elseif c_edge == 3
                if c_int == 0 return "p_on"
                elseif c_int == 4 count += 1 end
            elseif c_edge == 15
                x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
                if x_int > x count += 1
                elseif x_int == x return "p_on" end
            elseif (c_edge == 13) & ((c1==4) | (c2==4))
                    crossingTest(1,2,status,count)
            elseif (c_edge == 14) & ((c1==4) | (c2==4))
                    crossingTest(2,1,status,count)
            elseif c_edge == 7 count += 1
            elseif c_edge == 11 count = count
            elseif c_edge == 1
                if c_int == 0 return "p_on"
                elseif c_int == 4 crossingTest(1,2,status,count) end
            elseif c_edge == 2
                if c_int == 0 return "p_on"
                elseif c_int == 4 crossingTest(2,1,status,count) end
            elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
            elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
            elseif c_edge == 5
                if (c1==0) | (c2==0) return "p_on"
                else crossingTest(1,2,status,count) end
            elseif c_edge == 6
                if (c1==0) | (c2==0) return "p_on"
                else crossingTest(2,1,status,count) end
            elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
            elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
            end
        end
        if (round(count)%2)==1 
        	return "p_in"
        else 
        	return "p_out"
        end
    end
    return pointInPolygonClassification0
end


function crossRelation(FV::Array{Int64,2},EV::Array{Int64,2})
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = crossRelation(FW,EV)
end

function crossRelation(FV::Array{Array{Int64,1},1},EV::Array{Int64,2})
	sparseFV = cellComplex(FV)
	sparseEV = cellComplex(EV)
	sparseEF = (sparseFV' * sparseEV)'
	I,J,V = findnz(sparseEF)
	FE = [Int64[] for k=1:length(FV)]
	for (i,j,v) in zip(I,J,V) 
		if v==2
			push!(FE[j],i)
		end
	end
	return FE
end



function subModel(V::Array{Float64,2}, FV::Array{Int64,2}, EV::Array{Int64,2}, 
				bucket::Array{Int64,1},f::Int64,FE::Array{Array{Int64,1},1})
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = subModel(V,FW,EV,bucket,f,FE)
end

function subModel(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
					EV::Array{Int64,2},bucket::Array{Int64,1},f::Int64, FE::Array{Array{Int64,1},1})
	FW = [FV[f] for f in bucket]
	EW = hcat(collect(Set([EV[:,e] for f in bucket for e in FE[f]]))...)

	oldVerts = Set(Int64[])
	for k=1:length(bucket)
		union!(oldVerts, Set(FW[k]))
	end	
	olds = collect(oldVerts)
	vdict = Dict(zip(olds, 1:length(olds)))
	fz = [[vdict[v] for v in FW[k]] for k=1:length(FW)]
	ez = hcat([[vdict[v] for v in EW[:,k]] for k=1:size(EW,2)]...)
	z = hcat([V[:,v] for v in olds]...)
	pivot = find(bucket .== f)[1]
	return z,fz,ez,pivot
end


function visualize(Z,FZ::Array{Int64,2},EZ,f)
	FW = [FZ[:,k] for k=1:size(FZ,2)]
	out = visualize(Z,FW,EZ,f)
end

function visualize(Z,FZ::Array{Array{Int64,1},1},EZ,f)
	bucket = lar2hpc(Z,FZ)
	bucketEdges = lar2hpc(Z,EZ)
	params = PyObject(pyeval("list([1.,0.,0.,0.1,  0.,1.,0.,0.1,  
				0.,0.,1.,0.1, 0.,0.,0.,0.1, 100.])"))
	glass = p.MATERIAL(params)
	pivot = p.COLOR(p.RED)(p.JOIN(lar2hpc(Z,[FZ[f]])))
	p.VIEW(p.STRUCT([glass(bucket),bucketEdges,pivot]))
end


function meetZero( Z, vpair )
	(u,v) = vpair
    testValue = Z[:,u][3] * Z[:,v][3]
    if testValue > abs(10^-6.)
        return false
    else 
    	return true
    end
end

function crossZero( Z, vpair )
	(u,v) = vpair
    testValue = sign(vcode(Z[:,u])[3]) * sign(vcode(Z[:,v])[3])
    if testValue == -1
        return true
    else 
    	return false
    end
end

function spacePartition(V::Array{Float64,2}, FV::Array{Int64,2}, EV::Array{Int64,2},
						debug=false)
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = spacePartition(V,FW,EV,debug)
end


function spacePartition(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
						EV::Array{Int64,2},debug=false)
	""" input: face index f; candidate incident faces F[f]; """
	boxes = lar2boxes(V,FV)
	buckets = boxBucketing(boxes)
	FE = crossRelation(FV,EV)
	
	for (f,F) in enumerate(buckets)
		@show (f,F)
		""" F[f] submodel extraction from (V,FV,EV) """
		Z,FZ,EZ,pivot = subModel(V,FV,EV,F,f,FE)
		if debug visualize(Z,FZ,EZ,pivot) end
		
		""" Computation of submanifold map M moving f to z=0 """
		M = submanifoldMapping(V,FV,f)
		
		""" Transformation of F(f) by M, giving (W,EW,FW) := M(F(f)) """
		Z1 = vcat(Z,ones((1,size(Z,2))))
		Y = M * Z1
		Z = Y[ 1:3, : ]
		if debug visualize(Z,FZ,EZ,pivot) end
		
		""" filtering of EW edges traversing z=0, """
		fe = crossRelation(FZ,EZ)
		edges4face = [[EZ[:,e] for e in face] for (k,face) in enumerate(fe) if k!=pivot]
		dangling = [[vpair for vpair in edges if meetZero( Z, vpair )] 
					for edges in edges4face]
		single = Set(vcat([p.AA(sort)(pairs) for pairs in dangling if pairs!=[]]...))
		vpairs = collect(single)
		crossing = [[vpair for vpair in edges if crossZero( Z, vpair )] 
					for edges in edges4face]
		singlecross = Set(vcat([p.AA(sort)(pairs) for pairs in crossing if pairs!=[]]...))
		vpairscross = collect(singlecross)
		
		""" select pivot edges and 1D Lar model """
		pivotEdges = fe[pivot]
		v,ev = Z,[EZ[:,k] for k in pivotEdges]
		@show v
		@show ev

		""" Remove non-pivot vpairs with both external vertices """
		classify = pointInPolygonClassification(v[1:2,:],ev)
		function out(vk) 
			println("\neccomi")
			@show Z[1:2,vk], vk
			
			flag = classify(Z[1:2,vk])
			println(classify(Z[1:2,vk]))
			return flag != "p_out"
		end
		edges = [[v1,v2] for (v1,v2) in vpairs if (out(v1) | out(v2))]
		edges = hcat(vcat(edges,vpairscross)...)
		if debug visualize(Z,FZ,edges,pivot) end

		""" for each face in FZEZ, computation of the aligned set of points p(z=0) """
		
		
		""" Remove external vertices """
		
		""" Apply the inverse submanifold transform """
		
		""" Accumulate the submodel parts """
		
	end
	
	""" return the **valid** `LAR` 2-skeleton model `(W,FW,EW)` """
end


V,FV,EV = deepcopy((X,FX,EX))
