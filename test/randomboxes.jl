# Construction of independent buckets of containment boxes

# Given as input a list `randomLineArray` of pairs of 2D points, the function
# `containment2DBoxes` returns, in the same order, the list of _containment boxes_ of the
# input lines. A _containment box_ of a geometric object of dimension $d$ is defined as the
# minimal $d$-cuboid, equioriented with the reference frame, that contains the object. For a
# 2D line it is given by the tuple $(x1,y1,x2,y2)$, where $(x1,y1)$ is the point of minimal
# coordinates, and $(x2,y2)$ is the point of maximal  coordinates.

using LAR
include("src/inters.jl")

view(randomLines(1000,0.1)...)

randomLineArray = randomLines(200,0.2);
lines = lar2hpc(randomLineArray...);
yellow = p.COLOR(p.YELLOW)(lines);
boxes = containment2DBoxes(randomLineArray...);
rects = [box2rect(box) for box in boxes];
polylines = p.AA(polyline)(rects);
cyan = p.COLOR(p.CYAN)(p.STRUCT(polylines));
p.VIEW(p.STRUCT([yellow,cyan]))
