
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
viewLarIndices(W,EW,FW,0.75)

boxes = lar2boxes(W,FW)
parts = boxBucketing(boxes)



function larFromLines(datafile) # lineArray
	V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
	len = length(datafile)
	EV = collect(reshape(1:(len÷2), 2,(len÷4)))
	W,EW = lines2lar((V,EV))
	chains = boundary(W,EW)
	operator = boundaryOp(EW,chains)
	FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]
	W,FW,EW
end

V,FV,EV = larFromLines(datafile)



function chainAreas(V::Array{Float64,2},EV::Array{Int64,2},chains::Array{Int64,2})
	FE = [chains[:,f] for f=1:size(chains,2)]
	return chainAreas(V,EV,FE)
end


""" Implementation using integr.jl """
function chainAreas(V::Array{Float64,2}, EV::Array{Int64,2}, 
				chains::Array{Array{Int64,1},1})
V = vcat(V,zeros(1,size(V,2)))
pivots = [EV[:,abs(chain[1])][1] for chain in chains]
out = zeros(length(pivots))
for k=1:length(chains)
	area = 0
	triangles = [[] for h=1:length(chains[k])]
	for h=1:length(chains[k])
		edge = chains[k][h]
		v1,v2 = EV[:,abs(edge)]
		if sign(edge) == -1
			v1,v2=v2,v1
		end
		triangles[h] = Int[pivots[k],v1,v2]
	end
	P = V,triangles
	out[k] = surface(P,true)
end
return out
end


chainAreas(V,EV,chains)
