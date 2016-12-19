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


lineArray = randomLines(990,.3)

view(lineArray...)
V,EV = lines2lar(lineArray)
view(V,EV)
viewexploded(V,EV)
#viewLarIndices(V,EV,0.025)
chains = boundary(V,EV)
operator = boundaryOp(EV,chains)
FV = [sort(collect(Set(vcat([EV[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(V,EV,FV,0.1)

