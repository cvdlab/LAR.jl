# Test data
# using LAR


# Rounding of vectors to a given number of significant digits
vcode(v::Verts,prec=10^5) = zeros!(map(round, v*prec)/prec)
vcode(v::Vector,prec=10^5) = zeros!(map(round, v*prec)/prec)

# random 2D point with given number of digits.
# A single 2D random point, codified in floating point format, and with a fixed (quite small)
# number of digits, is returned by the \texttt{rand_point2d()} function, with no input
# parameters.
rand_point2d(prec=10^5) = vcode(rand(2,1))
    
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

# Generation of containment boxes of line segments  
# Input: a LAR pair V,EV
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
function centroid(boxes,xy)
	average = mean(boxes)
	n = Int(length(average)/2)
	median = (average[xy] + average[xy+n])/2
end

# Splitting the input above and below a median threshold
function splitOnThreshold(boxes,subset,coord)
    theBoxes = [boxes[k] for k in subset]
    threshold = centroid(theBoxes,coord)
    ncoords = length(boxes[1])
    a = coord % ncoords
    b = Int(a + ncoords/2)
    below,above = Int[],Int[]
    for k in subset
        if boxes[k][a] <= threshold 
        	push!(below, k) end 
    end
    for k in subset
        if boxes[k][b] >= threshold
        	push!(above, k) end 
    end
    Set{Int64}(below), Set{Int64}(above)
end

# Iterative splitting of a box array
function boxBuckets(boxes)
    bucket = Set(1:length(boxes))
    splittingStack = [bucket]
    finalBuckets = []
    
    # local function
	function splitting(bucket,below,above, finalBuckets,splittingStack)
		a = length(below)<4 & length(above) < 4
		b = length(setdiff(bucket,below)) < 4 
		c = length(setdiff(bucket,above)) < 4
		if ( a | b | c )
			push!(finalBuckets, [below]) 
			push!(finalBuckets, [above]) 
		else
			push!(splittingStack, below) 
			push!(splittingStack, above) 
		end
	end

    while splittingStack != []
        bucket = pop!(splittingStack)
        below,above = splitOnThreshold(boxes,bucket,1)
        below1,above1 = splitOnThreshold(boxes,above,2)
        below2,above2 = splitOnThreshold(boxes,below,2)
        splitting(above,below1,above1, finalBuckets,splittingStack)
        splitting(below,below2,above2, finalBuckets,splittingStack)      
    end
    finalBuckets = vcat(finalBuckets...)
end


