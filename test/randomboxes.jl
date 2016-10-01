# Construction of independent buckets of containment boxes

# Given as input a list `randomLineArray` of pairs of 2D points, the function
# `containment2DBoxes` returns, in the same order, the list of _containment boxes_ of the
# input lines. A _containment box_ of a geometric object of dimension $d$ is defined as the
# minimal $d$-cuboid, equioriented with the reference frame, that contains the object. For a
# 2D line it is given by the tuple $(x1,y1,x2,y2)$, where $(x1,y1)$ is the point of minimal
# coordinates, and $(x2,y2)$ is the point of maximal  coordinates.

using LAR

view(randomLines(10000,0.1)...)

randomLineArray = randomLines(200,0.2)
lines = p.AA(polyline)(linesFromLineArray(randomLineArray...))
boxes = containment2DBoxes(randomLineArray...)
rects = p.AA(box2rect)(boxes)
polylines = p.AA(polyline)(rects)
yellow = p.COLOR(p.YELLOW)(p.STRUCT(lines))
cyan = p.COLOR(p.CYAN)(p.STRUCT(polylines))
p.VIEW(p.STRUCT([yellow,cyan]))
