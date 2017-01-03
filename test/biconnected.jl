#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#

using LAR
include("src/inters.jl")
include("test/boundary.jl")

datafile = readcsv("test/svg/test2.lines");
#datafile = readcsv("test/svg/provaWall.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
len = length(datafile)
EV = collect(reshape(1:(len÷2), 2,(len÷4)))
view(V,EV)
lineArray = (V,EV)
W,EW = lines2lar(lineArray)
viewexploded(W,EW)
VV = vertices2vertices(W,EW)
viewLarIndices(W,EW,0.5)
#viewLarIndices(W,EW,0.025)
chains = boundary(W,EW)
operator = boundaryOp(EW,chains)
#println(full(operator))
FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(W,EW,FW,0.5)

boxes = lar2boxes(W,FW)
parts = boxBuckets(boxes)


lineArray = randomLines(500,.3)

view(lineArray...)
V,EV = lines2lar(lineArray)
view(V,EV)
viewexploded(V,EV)
#viewLarIndices(V,EV,0.025)
chains = boundary(V,EV)
operator = boundaryOp(EV,chains)
FV = [sort(collect(Set(vcat([EV[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(V,EV,FV,0.1)

boxes = lar2boxes(V,FV)
parts = boxBuckets(boxes)

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

V,EV,FV = boxes2lar(boxes)
view(V,EV)



