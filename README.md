# LAR.jl
### Geometric and topological modeling with chain complexes in Julia.


**Installation**: `julia> Pkg.clone("git://github.com/cvdlab/LAR.jl.git")`


## Basic Usage

```julia
using LAR
```

Generate and visualize a grid of 3-cubes. Poke around with the mouse pressed 
to interact.
```julia
v,cv = p.larCuboids((10,10,10));
viewexploded(v,cv)
```

Generate and visualize a grid of 2-cubes. Poke as before
```julia
v,fv = p.larCuboids((10,10));
viewexploded(v,fv)
```

Define a visualize a 2-mesh of convex cells. Poke as before
```julia
jsonmodel = """
	{"V" : [[5,0],[7,1],[9,0],[13,2],[15,4],[17,8],[14,9],[13,10],
	[11,11],[9,10],[5,9],[7,9],[3,8],[0,6],[2,3],[2,1],[8,3],
	[10,2],[13,4],[14,6],[13,7],[12,10],[11,9],[9,7],[7,7],[4,7],
	[2,6],[3,5],[4,2],[6,3],[11,4],[12,6],[12,7],[10,5],[8,5],
	[7,6],[5,5]],
	"FV" : [[0,1,16,28,29],[0,15,28],[1,2,17],[1,16,17,33],[2,3,17],
	[3,4,18,19],[3,17,18,30],[4,5,19],[5,6,19],[6,7,20,21,22,32],
	[6,19,20],[7,8,21],[8,9,21,22],[9,11,23,24],[9,22,23],
	[10,11,24,25],[10,12,25],[12,13,25,26],[13,14,27],[13,26,27],
	[14,15,28],[14,27,28,29,36],[16,29,34],[16,33,34],[17,30,33],
	[18,19,31],[18,30,31],[19,20,31,32],[22,23,32,33],[23,24,34,35],
	[23,33,34],[24,25,27,36],[24,35,36],[25,26,27],[29,34,35],
	[29,35,36],[30,31,32,33]]}
""";
model = json2larmodel(jsonmodel)
viewexploded(model.Verts',(model.Lar.FV))
```

```julia
lardict = JSON.parse(jsonmodel)
V = lardict["V"]
FV = lardict["FV"]
v,fv = p.larFacets((V,FV),2)
viewexploded(v,fv)
```


## Documentation

Very first basic implementation and interface with Python's `larlib`