using IntervalTrees
using LAR
include("src/inters.jl")
include("test/boundary.jl")

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


lineArray = randomLines(600,.3)

larview(lineArray...)
V,EV = lines2lar(lineArray)
larview(V,EV)
chains = boundary(V,EV)
operator = boundaryOp(EV,chains)
FV = [sort(collect(Set(vcat([EV[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(V,EV,FV,0.1)

boxes = lar2boxes(V,FV)
W,EW,FW = boxes2lar(boxes)
larview(W,EW)


buckets = boxBucketing(boxes)


cell = 1200
Z,EZ,_ = boxes2lar(lar2boxes(V,[FV[f] for f in buckets[cell]]))
pivot = p.COLOR(p.RED)(lar2hpc(boxes2lar(boxes[:,cell])[1:2]...))
bucket = lar2hpc(Z,EZ)
p.VIEW(p.STRUCT([bucket,pivot]))







