module mynewcode
function viewLarIndices(W,EW)
	submodel = lar2hpc(W,hcat(EW...))
	VV = PyObject(Any[Any[k] for k in  0:size(W,2)-1])
	EZ = map(Array{Int32},EW-1)
	EV = PyObject(Any[PyObject(EZ[e])[:tolist]() for e=1:length(EZ) ])
	V = PyObject(Any[PyObject(W[:,v])[:tolist]() for v=1:size(W,2) ])
	size = maximum(p.SIZE(Any[1,2])(submodel))/6
	hpc = p.larModelNumbering(1,1,1)(V,PyObject([VV,EV]),submodel,size)
	p.VIEW(hpc)
end
end
