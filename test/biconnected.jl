#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#

using LAR
include("src/inters.jl")

datafile = readcsv("test/svg/test1.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
len = length(datafile)
EV = collect(reshape(1:(len÷2), 2,(len÷4)))
view(V,EV)
lineArray = (V,EV)
W,EW = lines2lar(lineArray)
viewexploded(W,EW)
VV = vertices2vertices(W,EW)
viewLarIndices(W,EW,0.75)


lineArray = randomLines(300,.3)
view(lineArray...)
V,EV = lines2lar(lineArray)
view(V,EV)
viewexploded(V,EV)
viewLarIndices(V,EV,0.075)

















