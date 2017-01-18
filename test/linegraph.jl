using LAR
#include("src/inters.jl")

lineArray = randomLines(4000,0.1);
V,EV = lines2lar(lineArray);
viewexploded(V,EV)

lar = Lar();
lar.EV = [EV[:,k] for k=1:size(EV,2)];
ev = cellComplex(lar.EV);
ee = ev*ev';
EE = hcat(collect(findnz(ee))...)';
EW = Set([sort(EE[1:2,k]) for k=1:size(EE,2) if EE[3,k]==1 & (EE[1,k] != EE[2,k])]);
EW = hcat(collect(EW)...);
W = [mean([V[:,EV[1,h]] V[:,EV[2,h]]],2) for h=1:size(EV,2)];
W = hcat(W...);

primal = p.COLOR(p.CYAN)(lar2hpc(V,EV));
dual = p.COLOR(p.YELLOW)(lar2hpc(W,EW));
p.VIEW(p.STRUCT([primal,dual]))
