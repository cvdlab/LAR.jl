# Intersection of two line segments 
function segmentIntersect(boxes,lineArray,lineStorage)
	V,EV = lineArray
    function segmentIntersect0(h)
    	h1,h2 = EV[:,h]
        p1,p2 = V[:,h1],V[:,h2]
        line1 = vcode([p1 p2])
        (x1,y1),(x2,y2) = p1,p2
        B1,B2,B3,B4 = boxes[h]
        function segmentIntersect1(k)
			k1,k2 = EV[:,k]
			p3,p4 = V[:,k1],V[:,k2]
			line2 = vcode([p3 p4])
			(x3,y3),(x4,y4) = p3,p4
			b1,b2,b3,b4 = boxes[k]
            if ! |(b3<B1 , B3<b1 , b4<B2 , B4<b2)
                m23 = [p2 p3]
                m14 = [p1 p4]
                m = m23 - m14
                v = p3-p1
                a=m[1,1]; b=m[1,2]; c=m[2,1]; d=m[2,2]
                det = a*d-b*c
                if det != 0
                    m_inv = inv(m)
                    alpha, beta = m_inv*v
                    if (-0.0<=alpha<=1.0) & (-0.0<=beta<=1.0)
                        push!(lineStorage[line1], alpha)
                        push!(lineStorage[line2], beta)
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
        point = intersect1(line)
        if point != nothing
            push!(intersectionPoints, vcode(point))
        end
    end
    intersectionPoints
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
	boxes = containment2DBoxes(lineArray...)
	buckets = boxBuckets(boxes)
	intersectionPoints = Set{Array{Float64,1}}()
    for (h,bucket) in enumerate(buckets)
        pointBucket = lineBucketIntersect(boxes,lineArray, h,bucket, lineStorage)
        intersectionPoints = union(intersectionPoints,Set{Array{Float64,1}}(pointBucket))
    end
    frags = keys(lineStorage)
    params = [sort([x for x in Set(vcode(group))]) for group in values(lineStorage)]
    intersections = hcat([p for p in intersectionPoints]...)
    intersections,Dict(zip(frags,params))
end


include("src/inters.jl")
lineArray = randomLines(400,0.2)
V,EV = lineArray
[(V[:,EV[1,k]], V[:,EV[2,k]]) for k=1:size(EV,2)]
intPoints,lineFrags = lineIntersection(lineArray)


V = [[1.,0] [1.,1] [-1,0] [2,0] [0,-1] [1,-1]]
EV = [[1,5] [2,6] [3,4]]
lineArray = V,EV
intPoints,lineFrags = lineIntersection(lineArray)
