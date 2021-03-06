# Generation of lar model (V,EV) of a random set of 2D lines 

using LAR


V = [[0.,1] [1.,1] [-1,0] [2,0] [0,-1] [1,-1]]
EV = [[1,5] [2,6] [3,4]]
lineArray = V,EV
lineFrags = lineIntersection(lineArray)
W,EW = lines2lar(lineArray)
larview(W,EW)
viewexploded(W,EW)



lineArray = randomLines(200,0.4);
V,EV = lines2lar(lineArray);
hpcLines = lar2hpc(V,EV);
marker = p.CIRCLE(.002)([3,1]);
verts = PyObject(PyObject(V')[:tolist]());
markers = p.STRUCT(p.CONS(p.AA(p.T([1,2]))(verts))(marker));
redMarkers = p.COLOR(p.RED)(markers);
p.VIEW(p.STRUCT([hpcLines,redMarkers]))
viewexploded(V,EV)



