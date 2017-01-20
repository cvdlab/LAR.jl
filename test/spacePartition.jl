using LAR

#""" The `spacePartition` function takes as input a **non-valid** (with the meaning used in solid modeling field --- see~\cite{Requicha:1980:RRS:356827.356833}) `LAR` model of dimension d-1, i.e.~a triple `(V,FV,EV)`, The `spacePartition` function returns the **valid** `LAR` boundary model `(W,FW,EW)` of the space partition induced by `FV`. First an array `buckets` indexed on faces is computed. It contains the subsets of faces with greatest probability of intersecting each indexing face, respectively. 
#"""

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


function crossRelation(EV::Array{Int64,2},FV::Array{Array{Int64,1},1})
	sparseEV = cellComplex(EV)
	sparseFV = cellComplex(FV)
	sparseEF = (sparseFV' * sparseEV)'
	I,J,V = findnz(sparseEF)
	EF = [Int64[] for k=1:size(EV,2)]
	for (i,j,v) in zip(I,J,V) 
		if v==2
			push!(EF[i],j)
		end
	end
	return EF
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
					EV::Array{Int64,2},bucket::Array{Int64,1},f::Int64, 
					FE::Array{Array{Int64,1},1})
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

function edgefilter(Z,EZ,e)
	out = crossZero(Z,EZ[:,e]) | ((abs(Z[3,EZ[1,e]])<10^-6.) $ (abs(Z[3,EZ[2,e]])<10^-6.)) 
	return out
end

function intersectSegmentWithZero(p1::Array{Float64,1}, p2::Array{Float64,1})
	x1,y1,z1 = p1
	x2,y2,z2 = p2
	# z = z1 + alpha(z2 - z1)
	alpha_0 = (0 - z1)/(z2 - z1)
	point = Float64[x1 + alpha_0*(x2 - x1) , y1 + alpha_0*(y2 - y1)]
	return point
end


function larFromLines(datafile) # or lineArray
	V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
	len = length(datafile)
	EV = collect(reshape(1:(len÷2), 2,(len÷4)))
	W,EW = lines2lar((V,EV))
	chains = boundary(W,EW)
	operator = boundaryOp(EW,chains)
	FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]
	W,FW,EW
end


function remap(W,FW,EW)
	W,close,clusters,vmap = p.pruneVertices(W')	
	if typeof(W) == Array{Any,1} 
		W = Array{Float64,2}(hcat(W...))
	else 
		W = Array{Float64,2}(W')
	end
	vmap = vmap + 1
	EW = [sort([vmap[EW[1,k]],vmap[EW[2,k]]]) for k in 1:size(EW,2)]
	EW = collect(Set([EW[k] for k=1:length(EW)]))
	EW = hcat(EW...)
	FW = [sort([vmap[FW[k][h]] for h=1:length(FW[k])]) for k=1:length(FW)]
	FW = collect(Set([FW[k] for k=1:length(FW)]))
	return W,FW,EW
end


function checkLines4equalPoints(lines)
	out = []
	for line in lines
		if length(line) > 4
			b = reshape(line,2,Int(length(line)/2))
			line = vcat(p.pruneVertices(b')[1]...)
		end
		push!(out,line)
	end
	out
end


function spacePartition(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
						EV::Array{Int64,2},debug=false)
	""" input: face index f; candidate incident faces F[f]; """
	boxes = lar2boxes(V,FV)
	buckets = boxBucketing(boxes)
	FE = crossRelation(FV,EV)

	Vertices = Array{Float64,2}()
	Edges = Array{Int64,2}()
	Faces = Array{Array{Int64,1},1}()
	nverts = 0
	
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
		edges4face = [[e for e in face] for (k,face) in enumerate(fe) if k!=pivot]
		
		pivotEdges = fe[pivot]
		crossingEdges = [[e for e in edges if edgefilter(Z,EZ,e)] for edges in edges4face]
		edges = [edgeList for edgeList in crossingEdges if edgeList!=[]]
		EW = vcat(edges...,pivotEdges...)
		ez = hcat([EZ[:,e] for e in EW]...)
		if debug visualize(Z,FZ,ez,pivot) end

		""" for each face in EW, computate the aligned set of points p(z=0) """
		lines = Array{Array{Float64,1},1}()
		for face in edges
			#assert(length(face)==2) # may fail with non convex faces TODO: make general
			line = Array{Float64,1}()
			for e in face
				v1,v2 = EZ[:,e]
				px,py = vcode(intersectSegmentWithZero(Z[:,v1],Z[:,v2]))[1:2]
				push!(line,px); push!(line,py)
			end
			push!(lines,line)
		end
		for e in pivotEdges
			line = Array{Float64,1}()
			push!(line,vcode(Z[:,EZ[1,e]])[1]); push!(line,vcode(Z[:,EZ[1,e]])[2])
			push!(line,vcode(Z[:,EZ[2,e]])[1]); push!(line,vcode(Z[:,EZ[2,e]])[2])
			push!(lines,line)
		end
		lines = checkLines4equalPoints(lines)
		lineArray = hcat(lines...)
		W,FW,EW = larFromLines(lineArray')
		
		""" remove external cycle """
		chains = boundary(W,EW)
		if length(chains) == 2
			FW = [FW[1]]
		else
			areaPairing = chainAreas(W,EW,chains)
			maxAreas = maximum([abs(area) for area in areaPairing])
			FW = [face for (k,face) in enumerate(FW) if abs(areaPairing[k]) != maxAreas]
		end
		if debug viewLarIndices(W,EW,FW,3) end
				
		""" Apply the inverse submanifold transform """
		W1 = vcat(W,zeros((1,size(W,2))),ones((1,size(W,2))))
		Y = inv(M) * W1
		Z = Y[ 1:3, : ]
		if debug viewLarIndices(Z,EW,FW,3) end
		
		""" Accumulate the submodel parts """
		n = size(Z,2)
		if f==1
			Vertices = copy(Z)
			Edges = copy(EW)
			Faces = copy(FW)
			nverts = n
		else
			Vertices = hcat(Vertices,copy(Z))
			Edges = hcat(Edges,copy(EW)+nverts)
			Faces = vcat(Faces,copy(FW)+nverts)
			nverts += n
		end
	end
	""" return the **valid** `LAR` 2-skeleton model `(W,FW,EW)` """
	U,FU,EU = remap(Vertices,Faces,Edges)
end
















