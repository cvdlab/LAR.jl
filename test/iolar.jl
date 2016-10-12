using LAR
include("src/inters.jl")

lineArray = randomLines(400,0.3);
W,EW = lines2lar(lineArray);
model = W,EW

printLar("test/csv/test1", model)
V,EV = readLar("test/csv/test1")
view(V,EV)
viewexploded(V,EV)
