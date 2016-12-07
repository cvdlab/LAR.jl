#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#

using LAR
include("src/inters.jl")

function cols2any(EW)
	EZ = Any[]
	[push!(EZ, EW[:,k]) for k=1:size(EW,2)]
	EZ
end

datafile = readcsv("test/svg/test1.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2));
len = length(datafile);
EV = collect(reshape(1:(len÷2), 2,(len÷4)));
#view(V,EV);
lineArray = (V,EV)
W,EW = lines2lar(lineArray)
viewexploded(W,EW);
viewLarIndices(W,cols2any(EW),0.75)



lineArray = randomLines(300,.3)
Z,EZ = lines2lar(lineArray)
viewexploded(Z,EZ)
viewLarIndices(Z,cols2any(EZ),0.075)




viewexploded(V,EV);
VV = vertices2vertices(V,EV);
V,EW = biconnectedComponents(V,EV);
viewexploded(V,EW)
viewLarIndices(V,cols2any(EW),0.075)












