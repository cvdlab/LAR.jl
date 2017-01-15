
using LAR

datafile = readcsv("test/svg/test2.lines");
#datafile = readcsv("test/svg/provaWall.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
len = length(datafile)
EV = collect(reshape(1:(len÷2), 2,(len÷4)))
larview(V,EV)
lineArray = (V,EV)
W,EW = lines2lar(lineArray)
viewexploded(W,EW)
#VV = vertices2vertices(W,EW)
viewLarIndices(W,EW,0.75)
#viewLarIndices(W,EW,0.025)
chains = boundary(W,EW)
operator = boundaryOp(EW,chains)
#println(full(operator))
FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]
viewLarIndices(W,EW,FW,0.75)

boxes = lar2boxes(W,FW)
parts = boxBucketing(boxes)



V,FV,EV = larFromLines(datafile)



chainAreas(V,EV,chains)
