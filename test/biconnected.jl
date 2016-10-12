#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#

using LAR
include("src/inters.jl")

datafile = readcsv("test/svg/test1.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2));
len = length(datafile);
EV = collect(reshape(1:(len÷2), 2,(len÷4)));
view(V,EV);
lineArray = (V,EV);
W,EW = lines2lar(lineArray);
lar = Lar();
lar.EV = [EW[:,k] for k=1:size(EW,2)];
W,newcells,oldindex = larvalidate(W,lar);
EW = newcells.EV;
viewexploded(W,EW);


lineArray = randomLines(300,.3);
V,EV = lines2lar(lineArray);
viewexploded(V,EV);
VV = vertices2vertices(V,EV);
W,EW = biconnectedComponents(V,EV);
viewexploded(W,EW)

