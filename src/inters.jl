
# Rounding of vectors to a given number of significant digits
vcode(v::Verts,prec=10^5) = zeros!(map(round, v*prec)/prec)
vcode(v::Vector,prec=10^5) = zeros!(map(round, v*prec)/prec)

# random 2D point with given number of digits.
# A single 2D random point, codified in floating point format, and with a fixed (quite small)
# number of digits, is returned by the \texttt{rand_point2d()} function, with no input
# parameters.
rand_point2d(prec=10^15) = vcode(rand(2,1))
    
# Generation of a random line segment.
# A single random segment, scaled about its centroid by the \texttt{scaling} parameter, is
# returned by the \texttt{redge()} function, as a tuple ot two random points in the unit
# square.
function rand_edge(scaling=1)
	v1,v2 = (rand_point2d(), rand_point2d())
    c = (v1+v2)/2
    pos = rand_point2d()
    v1 = (v1-c)*scaling + pos
    v2 = (v2-c)*scaling + pos
    (vcode(v1), vcode(v2))
end

# Generation of random lines
# Return a LAR pair (V,EV)
function randomLines(lineNumber,scaling=1)
    L = hcat([hcat(rand_edge(scaling)...)  for k=1:lineNumber]...)
    xmin, ymin = min(L[1,:]...),min(L[2,:]...)
    v = [-xmin,-ymin]
    lines = [L[:,k]+v  for k=1:size(L,2)]
	V = hcat(lines...) 
	EV = reshape(collect(1:2*lineNumber),(2,lineNumber)) 
	V,EV
end

# Generation of containment boxes of LAR cells  
# Input: a LAR pair V,EV (or V,FV etc)
# Output: 
# Containment boxes 
function containment2DBoxes(V,EV)
	boxes = Vector{Float64}[]
	for k=1:size(EV,2) 
		(v1,v2) = EV[:,k]
		(x1,y1),(x2,y2) = V[:,v1],V[:,v2]
		box = [min(x1,x2),min(y1,y2),max(x1,x2),max(y1,y2)]
		boxes = push!(boxes,vcode(box))
	end
	boxes
end

function containment2DBoxes(V,EV)
	boxes = Vector{Float64}[]
	for k=1:size(EV,2) 
		(v1,v2) = EV[:,k]
		(x1,y1),(x2,y2) = V[:,v1],V[:,v2]
		box = [min(x1,x2),min(y1,y2),max(x1,x2),max(y1,y2)]
		boxes = push!(boxes,vcode(box))
	end
	boxes
end

# generate an HPC value from data of `Verts` type
polyline(verts::Verts) = p.POLYLINE(PyObject(map(PyObject,verts')))

# produce the 5 vertices of a closed rectangular polyline
# input is a box, i.e. a quadruple xmin,ymin,xmax,ymax()
function box2rect(box)
    x1,y1,x2,y2 = box
    verts = [[x1,y1] [x2,y1] [x2,y2] [x1,y2] [x1,y1]]
end

# linesFromLineArray(V,EV) = Any[[V[:,EV[:,k][1]] V[:,EV[:,k][2]]  ] for k=1:size(EV,2)]
# Implemented indirect indexing via multidimensional array
function linesFromLineArray(V,EV)
	result = Array(Float64, size(V,1), 2, size(EV,2))
	for k = 1 : size(EV,2)
		for j = 1 : 2
			for i = 1 : size(V,1)
				result[i,j,k] = V[i,EV[j,k]]
			end
		end
	end
	result
end

# Computation of the 1D centroid of a list of 2D boxes. The direction  
# is chosen depending on the value of the `xy in Set([1,2])` parameter. 
function centroid(boxes::Array{Float64,2},xyz)
	m = size(boxes,1)
	n = Int(m/2)
	average = zeros(m)
	for h=1:m
		average[h] = mean(boxes[h,:])
	end
	median = (average[xyz] + average[xyz+n])/2
end

function centroid(boxes::Array{Array{Float64,1},1},xyz)
	myboxes = hcat(boxes...)
	return centroid(myboxes::Array{Float64,2},xyz)
end


# Splitting the input above and below a median threshold
function splitOnThreshold(boxes,subset,coord)
	if subset==Set{Int64}()
		return Set{Int64}(), Set{Int64}()
	end
    theBoxes = [boxes[:,k] for k in subset]
    threshold = centroid(theBoxes,coord)
    ncoords = Int(floor(length(boxes[:,1])/2))
    a = coord % ncoords +1
    b = Int(a + ncoords)
    below,above = Int[],Int[]
    for k in subset
        if boxes[:,k][a] <= threshold 
        	push!(below, k) end 
    end
    for k in subset
        if boxes[:,k][b] >= threshold
        	push!(above, k) end 
    end
    Set{Int64}(below), Set{Int64}(above)
end

function geomPartitionate(boxes,buckets)
    geomInters = [Set{Int64}() for h=1:length(boxes)]
    for bucket in buckets
        for k in bucket
            geomInters[k] = union(geomInters[k],bucket)
        end
    end
    for (h,inters) in enumerate(geomInters)
        geomInters[h] = setdiff(geomInters[h],Set([h]))
    end
    geomInters
end

# local function
function splitting(bucket,below,above, finalBuckets,splittingStack)
	a = (length(below)<4) & (length(above)<4)
	b = length(setdiff(bucket,below)) < 7 
	c = length(setdiff(bucket,above)) < 7
	if ( a | b | c )
		push!(finalBuckets, below) 
		push!(finalBuckets, above) 
		println("pushed finalBuckets")
	else
		push!(splittingStack, below) 
		push!(splittingStack, above) 
		println("pushed splittingStack")
	end
end


function lar2boxes(V,CV::Array{Int64,2})
	m = size(V,1)*2
	boxes = zeros(Float64,(m,size(CV,2)))
	for k=1:size(CV,2)
	cell = CV[:,k]
		verts = hcat([V[:,v] for v in cell]...)
		pmin,pmax = Float64[],Float64[],Float64[]
		for h=1:size(verts,1)
			coords = verts[h,:]
			push!(pmin,minimum(coords))
			push!(pmax,maximum(coords))
		end
		box = vcat(pmin,pmax)
		boxes[:,k] = box
	end
	return boxes
end

function lar2boxes(V,CV::Array{Array{Int64,1},1})
	m = size(V,1)*2
	n = length(CV)
	boxes = zeros(Float64,(m,n))
	for k=1:n
		cell = CV[k]
		verts = hcat([V[:,v] for v in cell]...)
		pmin,pmax = Float64[],Float64[]
		for h=1:size(verts,1)
			coords = verts[h,:]
			push!(pmin,minimum(coords))
			push!(pmax,maximum(coords))
		end
		box = vcat(pmin,pmax)
		boxes[:,k] = box
	end
	return boxes
end



# Intersection of two line segments 
function segmentIntersect(boxes,lineArray,lineStorage)
	V,EV = lineArray
    function segmentIntersect0(h)
    	h1,h2 = EV[:,h]
        p1,p2 = V[:,h1],V[:,h2]
        line1 = vcode([p1 p2])
        #line1 = [p1 p2]
        (x1,y1),(x2,y2) = p1,p2
        B1,B2,B3,B4 = boxes[:,h]
        function segmentIntersect1(k)
			k1,k2 = EV[:,k]
			p3,p4 = V[:,k1],V[:,k2]
			line2 = vcode([p3 p4])
			(x3,y3),(x4,y4) = p3,p4
			b1,b2,b3,b4 = boxes[:,k]
            if ! |(b3<B1 , B3<b1 , b4<B2 , B4<b2)
                m23 = [p2 p3]
                m14 = [p1 p4]
                m = m23 - m14
                v = p3-p1
                a=m[1,1]; b=m[1,2]; c=m[2,1]; d=m[2,2]
                det = a*d-b*c
                if abs(det) > 10^-6. 
                    m_inv = inv(m)
                    alpha, beta = m_inv*v
                    if (-0.0<=alpha<=1.0) & (-0.0<=beta<=1.0)
                    	if 0.0 != vcode([alpha])[1] != 1.0
                        	push!(lineStorage[line1], alpha) 
                        end
                    	#if 0.0 != vcode([alpha])[1] != 1.0
                    	if 0.0 != vcode([beta])[1] != 1.0
	                        push!(lineStorage[line2], beta)
	                    end
                        return p1+alpha*(p2-p1)
                    end
                end
            end
            #return nothing
        end
    end
end



# Brute force intersection within the buckets
function lineBucketIntersect(boxes,lineArray, h,bucket, lineStorage)
    intersect0 = segmentIntersect(boxes,lineArray,lineStorage)
    intersectionPoints = Array{Float64,1}[]
    intersect1 = intersect0(h)
    for line in bucket
    	#@show line
        point = intersect1(line)
        if point != nothing
            push!(intersectionPoints, vcode(point))
        end
    end
    intersectionPoints
end

function boxBucketing(boxes::Array{Float64,2})
	nboxes = size(boxes,2)
	dim = Int(size(boxes,1)/2)
	trees,boxes1D = Any[],Any[]
	# preparation of d interval-trees
	for d=1:dim
		intervals = [ (boxes[d,k], boxes[d+dim,k], k) for k=1:nboxes ];
		push!(boxes1D, intervals)
		tree1D = IntervalTree{Float64, IntervalValue{Float64, Int64}}();
		for triple in intervals
			push!(tree1D, IntervalValue{Float64, Int64}(triple...));
		end
		push!(trees, tree1D)
	end
	# execution of spatial queries
	buckets,queries = Any[],Any[]
	for k=1:nboxes
		sets = Any[]
		for d=1:dim
			push!(sets, Set(Int64[]))
			query = boxes1D[d][k]
			for item in intersect(trees[d],query[1:2])
				push!(sets[d], item.value)
			end
		end
		push!(buckets,intersect(sets...))
	end
	buckets = [sort(collect(bucket)) for bucket in buckets]
	return buckets
end


# Accelerate intersection of lines 
function lineIntersection(lineArray)
	V,EV = lineArray
	lineStorage = Dict{Array{Float64,2},Array{Float64,1}}()
	for k = 1:size(EV,2)
		v1,v2 = EV[:,k]
		p1,p2 = V[:,v1],V[:,v2]
		key = vcode([p1 p2])
		lineStorage[key] = Int[]
	end
	boxes = lar2boxes(lineArray...)
	buckets = boxBucketing(boxes)
    for (h,bucket) in enumerate(buckets)
    	#@show h,bucket
        pointBucket = lineBucketIntersect(boxes,lineArray, h,bucket, lineStorage)
    end
    frags = keys(lineStorage)
    params = [sort([x for x in Set(vcode(group))]) for group in values(lineStorage)]
    lineFrags = Dict(zip(frags,params))
end

function check4DoubleEdges(EV)
	for k=1:(size(EV,2)-1)
		if EV[:,k]==EV[:,k+1]
			
		end
	end
end


function larModelCheck(V,EV)
	vertexIndex,k = Dict{Int64,Int64}(),0
	W = sort(collect(Set(EV)))
	for (k,v) in enumerate(W)
		vertexIndex[v] = k
	end
	EW = deepcopy(EV)
	for k=1:size(EV,2)
		EW[1,k] = vertexIndex[EV[1,k]]
		EW[2,k] = vertexIndex[EV[2,k]]
	end
	vkeys = p.TRANS(sort(collect(vertexIndex)))'[:,1]
	V = hcat([V[:,k] for k in vkeys]...)
	EV = hcat([sort(EW[:,k]) for k=1:size(EW,2)]...)
	EV = sortcols(EV, by=x->(x[1],x[2]))
	EV = sortrows(EV)
	EV = check4Edges(EV)
	V,EV
end

# Transform a lineArray (Array of pairs of 2D points) into a 1-complex (V,EV)
function lines2lar(lineArray,prec=10^15)
	lineFrags = lineIntersection(lineArray)
	
	verts = zeros(2,1)
	index = 0
	EV = Array{Any,1}()
	for k in 1:length(lineFrags)
		ps = collect(keys(lineFrags))[k]
		alphas = vcat([[0.],collect(values(lineFrags))[k],[1.]]...)
		pts = [ (ps[:,1]*(1.-alphas[i])+ps[:,2]*alphas[i]) for i=1:length(alphas) ]
		verts = hcat(verts,hcat(pts...))
		for h=(index+1):(index+length(pts)-1)
			push!(EV, [h,h+1])
		end
		index += length(pts)
	end
	V,EV = verts[:,2:end],hcat([EV[k] for k=1:length(EV)]...)

	W,close,clusters,vmap = p.pruneVertices(V')	
	if typeof(W) == Array{Any,1} 
		Z = Array{Float64,2}(hcat(W...))
	else 
		Z = Array{Float64,2}(W)
	end
	vmap = vmap + 1
	EW = [[vmap[EV[1,k]],vmap[EV[2,k]]] for k in 1:size(EV,2)]
	EZ = hcat(EW...)
	V,EV = biconnectedComponents(Z,EZ)
	if size(V,2)==2 V=V' end
	V,EV = larModelCheck(V,EV)
end



# Adjacency lists of 1-complex vertices 
function vertices2vertices(V,EV)
	lar = Lar()
	lar.EV = [EV[:,k] for k=1:size(EV,2)]
	ev = cellComplex(lar.EV)
	vv = ev'*ev
	triples = hcat(findnz(vv)...)'
	VV = [Int64[] for k=1:size(vv,1)]
	for k=1:size(triples,2)
		row,col,datum = triples[:,k]
		if row != col
			push!(VV[col],row)
		end
	end
	return VV
end


# Main procedure for biconnected components 
function biconnectedComponents(W,EV)
    V = 1:length(W)
    count = 0
    stack,out = [],[]
    visited = [false for v in V]
    parent = [-1 for v in V]
    d = [-1 for v in V]
    low = [-1 for v in V]
    VV = vertices2vertices(W,EV)
    for u in V
        if !visited[u]
            DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
        end
    end
    out = [component for component in out if length(component) > 1]
    EV = [[u,v] for (u,v) in vcat(out...)]
    W,hcat(EV...)
end



# Hopcroft-Tarjan algorithm 
function DFV_visit( VV,out,count,visited,parent,d,low,stack,u )
    visited[u] = true
    count += 1
    d[u] = count
    low[u] = d[u]
    if u > length(VV) return() end
    for v in VV[u] 
        if !visited[v]
            push!(stack, (u,v))   # possible bug: []
            #println("stack = $stack")
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]
                push!(out,outputComp(stack,u,v))
                #println("\tout = $out")
                #println("\tv = $v")
                if v == length(VV) return() end
            end
            low[u] = min( low[u], low[v] )
        else
            if ! (parent[u]==v) & (d[v] < d[u])
                push!(stack, (u,v))   # possible bug: []
                low[u] = min( low[u], d[v] )
            end
        end
    end
end

# Output of biconnected components 
function outputComp(stack,u,v)
    out = []
    while true
        elem = pop!(stack)
        push!(out, elem)
        if elem == (u,v) 
        	break
        end
    end
    out
end



function larFromLines(datafile) # lineArray'
	V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
	len = length(datafile)
	EV = collect(reshape(1:(len÷2), 2,(len÷4)))
	W,EW = lines2lar((V,EV))
	chains = boundary(W,EW)
	operator = boundaryOp(EW,chains)
	FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]
	W,FW,EW
end
