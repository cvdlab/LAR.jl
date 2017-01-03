using IntervalTrees

using LAR
include("src/inters.jl")
include("test/boundary.jl")

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
	return [sort(collect(bucket)) for bucket in buckets]
end


lineArray = randomLines(600,.3)

view(lineArray...)
V,EV = lines2lar(lineArray)
view(V,EV)
chains = boundary(V,EV)
operator = boundaryOp(EV,chains)
FV = [sort(collect(Set(vcat([EV[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(V,EV,FV,0.1)

boxes = lar2boxes(V,FV)

function boxes2lar(boxes)
	V = Array{Float64,1}[]
	EV = Array{Int64,1}[]
	FV = Array{Int64,1}[]
	for k=1:size(boxes,2)
		xm,ym,xM,yM = boxes[:,k]
		push!(V,[xm,ym])
		push!(V,[xM,yM])
		push!(V,[xm,yM])
		push!(V,[xM,ym])
		push!(EV,[ 4(k-1)+1, 4(k-1)+3 ])
		push!(EV,[ 4(k-1)+1, 4(k-1)+4 ])
		push!(EV,[ 4(k-1)+2, 4(k-1)+3 ])
		push!(EV,[ 4(k-1)+2, 4(k-1)+4 ])
		push!(FV,[ 4(k-1)+1, 4(k-1)+2, 4(k-1)+3, 4(k-1)+4])
	end
	return hcat(V...), hcat(EV...), hcat(FV...)
end

W,EW,FW = boxes2lar(boxes)
view(W,EW)


buckets = boxBucketing(boxes)


cell = 1200
Z,EZ,_ = boxes2lar(lar2boxes(V,[FV[f] for f in buckets[cell]]))
pivot = p.COLOR(p.RED)(lar2hpc(boxes2lar(boxes[:,cell])[1:2]...))
bucket = lar2hpc(Z,EZ)
p.VIEW(p.STRUCT([bucket,pivot]))







