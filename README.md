# LAR.jl
### Geometric and topological modeling with chain complexes in Julia.

**Pre-requisite**: first install [pyplasm](https://github.com/plasm-language/pyplasm/blob/master/README.rst) and [larlib](https://pypi.python.org/pypi/larlib/) for Python 2.7

**Installation**: `julia> Pkg.clone("git://github.com/cvdlab/LAR.jl.git")`


## Basic Usage

```julia
using LAR
```

Generate and visualize a grid of 3-cubes. Poke around with the mouse pressed 
to interact.
```julia
v,cv = p.larCuboids((10,10,10));
viewexploded(v',(cv+1)')
```

Generate and visualize a grid of 2-cubes. Poke as before
```julia
v,fv = p.larCuboids((10,10));
viewexploded(v',(fv+1)')
```

Define a visualize a 2-mesh of convex cells. Poke as before
```julia
jsonmodel = """
	{"V" : [[5.0,0.0],[7.0,1.0],[9.0,0.0],[13.0,2.0],[15.0,4.0],[17.0,
	8.0],[14.0,9.0],[13.0,10.0],[11.0,11.0],[9.0,10.0],[5.0,9.0],[7.0,
	9.0],[3.0,8.0],[0.0,6.0],[2.0,3.0],[2.0,1.0],[8.0,3.0],[10.0,2.0],
	[13.0,4.0],[14.0,6.0],[13.0,7.0],[12.0,10.0],[11.0,9.0],[9.0,7.0],
	[7.0,7.0],[4.0,7.0],[2.0,6.0],[3.0,5.0],[4.0,2.0],[6.0,3.0],[11.0,
	4.0],[12.0,6.0],[12.0,7.0],[10.0,5.0],[8.0,5.0],[7.0,6.0],[5.0,5.0]],
	"FV" : [[0,1,16,28,29],[0,15,28],[1,2,17],[1,16,17,33],[2,3,17],
	[3,4,18,19],[3,17,18,30],[4,5,19],[5,6,19],[6,7,20,21,22,32],
	[6,19,20],[7,8,21],[8,9,21,22],[9,11,23,24],[9,22,23],
	[10,11,24,25],[10,12,25],[12,13,25,26],[13,14,27],[13,26,27],
	[14,15,28],[14,27,28,29,36],[16,29,34],[16,33,34],[17,30,33],
	[18,19,31],[18,30,31],[19,20,31,32],[22,23,32,33],[23,24,34,35],
	[23,33,34],[24,25,27,36],[24,35,36],[25,26,27],[29,34,35],
	[29,35,36],[30,31,32,33],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]]}
""";
model = json2larmodel(jsonmodel);
viewexploded(model.Verts,rebase(model.Lar.FV[1:end-1]))
```

```julia
lardict = JSON.parse(jsonmodel);
V = lardict["V"];
FV = lardict["FV"];
v,ev = p.larFacets((V,FV),2);
viewexploded(v',(ev+1)')
ev = Any[ev[k,:]+1 for k=1:size(ev,1)]
viewLarIndices(v',ev)
```


## Documentation

Very first basic implementation and interface with Python's `larlib`