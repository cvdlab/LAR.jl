# module Lar-core
# basic implementation of lar-core.jl module

using JSON

using PyCall

@pyimport larlib as p

# root of LAR type hierarchy
abstract LarModel

# top LAR subtypes
abstract LarGeometry <: LarModel
abstract CellComplex <: LarModel
abstract ChainComplex <: LarModel

# concrete type for models
type Model <: LarModel
	Verts
	Lar
end

# `Cells` type is a list (of lists of integers)
typealias Cells Array{Any,1}

# concrete type `LAR` is a quadruple of `Cells`
type Lar <: CellComplex
    VV::Cells
    EV::Cells
    FV::Cells
    CV::Cells
end

# constructor of an empty quadruple of lists
function Lar() 
	Lar(Any[],Any[],Any[],Any[]) 
end

# constructor of a `Lar` instance (some fields possibly empty)
function Lar(bases::Array{Any,1}) 
	Lar( vv::Array{Any,2}, ev::Array{Any,2}, fv::Array{Any,2}, cv::Array{Any,2} ) 
end

# `Basis` type as alias for Sparse binay Matrix by Columns
typealias Basis SparseMatrixCSC

# `ChainBases` is a quadruple (just 3D, for now) of `Basis` binary matrices 
type ChainBases <: ChainComplex
    M0::Basis
    M1::Basis
    M2::Basis
    M3::Basis
end

# constructor of an empty `ChainBases`
function ChainBases() 
	m0 = cellComplex()
	m1 = cellComplex()
	m2 = cellComplex()
	m3 = cellComplex()
	ChainBases(m0,m1,m2,m3) 
end

# constructor of a Basis instance from its cells
function cellComplex(cells::Cells)
	I,J,V = Int[],Int[],Int[]
	for (i,row) in enumerate(cells)
		for (j,col) in enumerate(row)
			push!(I,i); push!(J,col); push!(V,1)
		end
	end
	sparse(I,J,V)
end

function cellComplex(cells::Array{Array{Int64,1},1})
	I,J,V = Int[],Int[],Int[]
	m,n = length(cells),maximum(vcat(cells...))
	for j=1:m
		for i=1:length(cells[j])
			push!(I,cells[j][i]); push!(J,j); push!(V,1)
		end
	end
	sparse(I,J,V)
end

# constructor of a Basis instance from its cells
function cellComplex(cells::Array{Any,2})
	I,J,V = Int[],Int[],Int[]
	m,n = size(cells)
	for i=1:m
		for j=1:n
			push!(I,cells[i,j]); push!(J,j); push!(V,1)
		end
	end
	sparse(I,J,V)
end

# constructor of a Basis instance from its cells
function cellComplex(cells::Array{Int64,2})
	I,J,V = Int[],Int[],Int[]
	m,n = size(cells)
	for i=1:m
		for j=1:n
			push!(I,cells[i,j]); push!(J,j); push!(V,1)
		end
	end
	sparse(I,J,V)
end

# constructor of an empty Basis instance
function cellComplex()
	sparse(Int[],Int[],Int[])
end

# alias type used for boundary and coboundary operators
typealias Map SparseMatrixCSC

# concrete type `CCMap` (for Chain-Cochain Map) for chain complexes (3D)
type CCMap <: ChainComplex
    B1::Map
    B2::Map
    B3::Map
end

# constructor of empty `CCMap`
function CCMap() 
	CCMap(sparse(Int[],Int[],Int[])) 
end

# alias `Verts` for array of vertices in any embedding space (nD)
typealias Verts Array{Float64,2}

# Constructor of vertices from a Lar dictionary
function Verts(lar::Dict)
	verts = get(lar, "V", 0)
	if verts != 0
		if length(size(verts))==1
			verts = map(Float64, reduce(hcat,verts))
			vs = Verts(verts) 
		else 
			println("\nerror: vertex input")
		end
	end
end

# conversion from a JSON dictionary to a concrete `Model` type,
# i.e. to a pair of (Verts, Lar) types
function json2larmodel(strng::AbstractString)
	lardict = JSON.parse(strng)
	vs = Verts(lardict)
	verts = get(lardict, "VV", 0)
	edges = get(lardict, "EV", 0)
	faces = get(lardict, "FV", 0)
	cells = get(lardict, "CV", 0)
	lar = Lar()
	if verts != 0 lar.VV = verts end
	if edges != 0 lar.EV = edges end
	if faces != 0 lar.FV = faces end
	if cells != 0 lar.CV = cells end
	Model(vs,lar)
end

# used to convert a 2-array to a list of lists
function listOfList(cellArray::Array{Any,2})
	out = []
	for k=1:size(EV1)[1]
		cell = copy(cellArray[k,:])
		cell = cell[cell .!= ""]
		push!(out, cell)
	end
	out
end

# transform a 0-based array to a 1-based array
function rebase(faces::Array{Any,1},shift=+1)
	Any[Int[v+shift for v in face] for face in faces]
end

# transform a 0-based array to a 1-based array
function rebase(faces::Array{Any,2},shift=+1)
	faces + shift
end

# transform a `larlib` pair (varts,bases) into a Julia triple
# of instances of (Verts,Lar,ChainComplex) types
function importmodel(larmodel)
	if typeof(larmodel) == Tuple{Array{Any,2},Array{Any,1}}
		verts,bases = larmodel
		n = length(verts)
		chainbases = ChainBases()
		cells = Lar()
		v =  map(Float64,verts)'
		if length(bases) in (3,4)
			vv = rebase(bases[1]')
			ev = rebase(bases[2]')
			fv = rebase(bases[3]') 
			cells.VV = Vector{Int}[vv[:,i] for i=1:size(vv,2)]
			cells.EV = Vector{Int}[ev[:,i] for i=1:size(ev,2)]
			cells.FV = Vector{Int}[fv[:,i] for i=1:size(fv,2)]
			chainbases.M0 = cellComplex(vv)
			chainbases.M1 = cellComplex(ev)
			chainbases.M2 = cellComplex(fv) 
		end
		if length(bases) == 4
			cv = rebase(larmodel[2][4]')
			cells.CV = Vector{Int}[cv[:,i] for i=1:size(cv,2)]
			chainbases.M3 = cellComplex(cv)
		end
		v,cells,chainbases
	end
end

function lar2hpc(V::Array{Float64,2}, EV::Array{Any,1})
	EV = map(Int32, EV)
	lar2hpc(V,EV)
end

function lar2hpc(V::Array{Float64,2}, EV::Array{Array{Int64,1},1})
	EV = hcat(EV...)
	lar2hpc(V,EV)
end

# function lar2hpc(V::Array{Float64,2}, EV::Array{Array{Int64,1},1})
# 	EV = map(Any, EV)
# 	lar2hpc(V,EV)
# end

function lar2hpc(V::Array{Float64,2}, EV::Array{Int32,2})
	a,b = PyObject(V'), PyObject(EV')
	verts = PyObject(a[:tolist]())
	cells = PyObject(b[:tolist]())
	p.MKPOL([verts,cells,1])
end

function lar2hpc(V::Array{Float64,2}, EV::Array{Int64,2})
	EV = map(Int32, EV)
	a,b = PyObject(V'), PyObject(EV')
	verts = PyObject(a[:tolist]())
	cells = PyObject(b[:tolist]())
	p.MKPOL([verts,cells,1])
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Int64,2})
	p.VIEW(lar2hpc(V,EV))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Any,2})
	EV = map(Int64, EV)
	view(V,EV)
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Any,2}, EV::Array{Any,1})
	a,b = PyObject(V'), PyObject(EV)
	p.VIEW(p.MKPOL([a,b,1]))
end


# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Any,1})
	p.VIEW(lar2hpc(V,EV))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Array{Int64,1},1})
	a,b = PyObject(V'), PyObject(EV)
	verts = PyObject(a[:tolist]())
	cells = b
	p.VIEW(p.MKPOL([verts,cells,1]))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Any,2}, EV::Array{Any,2})
	a,b = PyObject(V'), PyObject(EV')
	verts = a
	cells = b
	p.VIEW(p.MKPOL([verts,cells,1]))
end

# generate the triple of scaling factors needed by `larlib` explode
function scalingargs(scaleargs)
	if length(scaleargs)==1
		sx = sy = sz = scaleargs[1]
	elseif length(scaleargs)==3
		sx, sy, sz = scaleargs
	else
		println("error: wrong number of scaling arguments")
	end
	sx, sy, sz
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Float64,2}, fv::Array{Array{Int64,1},1}, scaleargs=1.2)
	sx, sy, sz = scalingargs(scaleargs)
	a = PyObject(v')
	b = PyObject(Array{Any,1}[fv[k]-1 for k in 1:length(fv)])
	verts = PyObject(a[:tolist]())
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((verts,b))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Float64,2}, fv::Array{Any,1}, scaleargs=1.2)
	sx, sy, sz = scalingargs(scaleargs)
	a = PyObject(v')
	b = PyObject(Array{Any,1}[fv[k]-1 for k in 1:length(fv)])
	verts = PyObject(a[:tolist]())
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((verts,b))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Any,2}, fv::Array{Any,1}, scaleargs=1.2)
	sx, sy, sz = scalingargs(scaleargs)
	verts = PyObject(v')
	if typeof(verts)==PyObject
		a = verts
	else
		a = PyObject(verts[:tolist]())
	end
	b = PyObject(Array{Any,1}[fv[k]-1 for k in 1:length(fv)])
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((a,b))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Int64,2}, fv::Array{Int64,2}, scaleargs=(1.2,))
	v = map(Float64,v)
	viewexploded(v, fv, scaleargs)
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Float64,2}, fv::Array{Int64,2}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	a,b = PyObject(v'), PyObject(fv'-1)
	verts = PyObject(a[:tolist]())
	cells = PyObject(b[:tolist]())
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((verts,cells))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Float64,2}, fv::Array{Any,2}, scaleargs=(1.2,))
	fv = map(Int64, fv)
	viewexploded(v, fv, scaleargs)
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Any,2}, fv::Array{Any,2}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	fv = map(Int,fv)-1
	cells = [(fv[:,k]) for k in 1:size(fv,2)]
	a,b = PyObject(v'), PyObject(cells)
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((a,b))))
end

function viewexploded(v::Array{Any,2}, fv::Array{Int64,2}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	fv = fv-1
	cells = [(fv[:,k]) for k in 1:size(fv,2)]
	a,b = PyObject(v'), PyObject(cells)
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((a,b))))
end

# embed vertices array in a higher dimensional space, adding `ndims` zero coords
function embed(verts::Verts; ndims=1)
	vs = copy(verts)
	thetype = typeof(vs[1,1])
	added = zeros(thetype, (ndims, size(vs)[2]))
	vs = [vs ; added]
	out = Verts(vs)
end

# translate the columns of `V` matrix by sum with `t` vector
function translate( t, V )
	broadcast(+,t,V)
end

# scale the columns of `V` matrix by product times `s` vector
function scale( s, V )
	broadcast(*,s,V)
end

# rotate the columns of `V` matrix by properly using the `args` parameters
function rotate(args,V)
	n = length(args)
	if n == 1 # rotation in 2D
		angle = args[1]; Cos = cos(angle); Sin = sin(angle)
		mat = eye(2)
		mat[1,1] = Cos;    mat[1,2] = -Sin;
		mat[2,1] = Sin;    mat[2,2] = Cos;
	elseif n == 3 # rotation in 3D
		mat = eye(3)
		args = PyObject(PyObject(args)[:tolist]())
		angle = p.VECTNORM(args); axis = p.UNITVECT(args)
		Cos = cos(angle); Sin = sin(angle)
		# elementary rotations (in 3D)
		if axis[2]==axis[3]==0.0    # rotation about x
			mat[2,2] = Cos;    mat[2,3] = -Sin;
			mat[3,2] = Sin;    mat[3,3] = Cos;
		elseif axis[1]==axis[3]==0.0    # rotation about y
			mat[1,1] = Cos;    mat[1,3] = Sin;
			mat[3,1] = -Sin;    mat[3,3] = Cos;
		elseif axis[1]==axis[2]==0.0    # rotation about z
			mat[1,1] = Cos;    mat[1,2] = -Sin;
			mat[2,1] = Sin;    mat[2,2] = Cos;
		# general rotations (in 3D)
		else  # general 3D rotation (Rodrigues' rotation formula)    
			I = eye(3) ; u = axis;
			Ux = Array([
				 0        -u[3]      u[2];
				 u[3]        0      -u[1];
				-u[2]      u[1]        0 ])
			UU = Array([
				 u[1]*u[1]    u[1]*u[2]    u[1]*u[3];
				 u[2]*u[1]    u[2]*u[2]    u[2]*u[3];
				 u[3]*u[1]    u[3]*u[2]    u[3]*u[3] ])
			mat = Cos*I + Sin*Ux + (1.0-Cos)*UU
		end
	end
	mat*V
end


# Identify -0.0 values with 0.0 values
function zeros!(W)
	for k=1:length(W)
	   if W[k] == -0.0 W[k] = 0.0 end
	end
	W
end


# vectorized parametric toroidal surface
function toroidalmap(V,r,R,m=40,n=80,angle1=2pi,angle2=2pi)
	V = map(Float64,V)
	V = (scale(V, [angle1/m angle2/n]))
	x =  (R .+ r.*cos(V[:,1])) .* cos(V[:,2])
	y =  (R .+ r.*cos(V[:,1])) .* sin(V[:,2])
	z = -r .* sin(V[:,1])
	[x y z]'
end


# valid (closed, i.e. no boundary) toroidal surface
function toroidalsurface(r,R,m=10,n=20,angle1=2pi,angle2=2pi)
	larmodel = p.larCuboids((m,n),true)
	verts,cells,chainbases = importmodel(larmodel)
	V = toroidalmap(verts',r,R,m,n)
	W,newcells,oldindex = larvalidate(V,cells)
	Z = map(Float64,hcat([V[:,k] for k in oldindex]...))
	Z,newcells
end

# topological validation of cellular complexes 
function larvalidate(V,lar,prec=10^5) # prec for precision
	W = zeros!(map(round, V*prec)/prec)
	vdict = Dict([ ("$(W[:,k])",k) for k=1:size(W,2) ])
	verts = Dict(zip( keys(vdict), 1:length(vdict) ))
	newindex = convertindex(W,verts,prec)
	oldindex = invertindex(newindex)
	vs = [eval(parse(key))' for key in keys(vdict)]
	W = vcat(vs...)'
	newcells = Lar()
	if lar.FV != []
		FW = relink(lar.FV,newindex)
		newcells.FV = [ FW[:,k][:] for k=1:size(FW,2) ] 
	end
	if lar.EV != []
		EW = relink(lar.EV,newindex)
		newcells.EV = [ EW[:,k][:] for k=1:size(EW,2) ]
	end
	if lar.VV != []
		VW = relink(lar.VV,newindex)
		newcells.VV = [ VW[:,k][:] for k=1:size(VW,2) ]
	end
	W,newcells,oldindex
end

# relinking of chain basis (d-cells, named FV) after topology validation
# including removal of multiple vertices (TODO ...)
# and 
function relink(basis,newindex)
	m0,n0 = length(basis),length(basis[1])
	FV = reshape(vcat(basis...),n0,m0)
	FW = round(Int64, zeros(size(FV)))
	for k in 1:size(FV,2)
		for h in 1:size(FV,1)
			FW[h,k] = newindex[FV[h,k]]
		end
		FW[:,k] = sort(FW[:,k])
	end
	FW
end

# mapping between vertex indexing in dictionary and old indexing
function convertindex(W,vdict,precision)
	[vdict[ string(zeros!(map(round,W[:,k]*precision)/precision)) ] for k=1:size(W,2)]
end

function invertindex(index)
	n = max(index...)
	out = round(Int,zeros(n))
	for k=1:length(index)
		out[index[k]] = k
	end
	out
end

# Printing a Lar model (V,EV)
function printLar(modelname, model)
	V,EV = model
	
	function writefile(filename,data::Matrix)
		out_file = open(filename, "w")
		n,m = size(data)
		# file records printed by column
		for h=1:m
			for k=1:n-1
				print(out_file, data[k,h], ", ")
			end
			println(out_file, data[end,h])
		end
		close(out_file)
	end
	
	# vertex file
	filename = join([modelname,".V"])
	writefile(filename,V)
	# cell file
	filename = join([modelname,".EV"])
	writefile(filename,EV)
end

# Reading a Lar model (V,EV)
function readLar(modelname)
	filename = join([modelname,".V"])
	W = readcsv(filename)
	filename = join([modelname,".EV"])
	EW = readcsv(filename,Int64)
	W',EW'
end

function cols2any(EW)
	EZ = Any[]
	[push!(EZ, EW[:,k]) for k=1:size(EW,2)]
	EZ
end

# Produce an indexed view of Lar model (W,EW)
function viewLarIndices(W,EW,unit=1.0)
	EW = cols2any(EW)
	function convertData(W,EW)
		# Shifted vertices
		Z = zeros(length(W)+size(W,1))
		for k=1:length(W)
			Z[k + size(W,1)] = W[k]
		end
		Z = reshape(Z,size(W,1),size(W,2)+1)
		# shifted edges
		EZ = Any[[1,1]]
		for k=1:length(EW)
			push!(EZ, EW[k]+1)
		end
		Z,EZ
	end
	V,EV = convertData(W,EW)
	submodel = lar2hpc(V,hcat(EV...))
	VV = PyObject(Any[Any[k] for k in  0:size(V,2)-1])
	EZ = map(Array{Int32},EV-1)
	EV = PyObject(Any[PyObject(EZ[e])[:tolist]() for e=1:length(EZ) ])
	V = PyObject(Any[PyObject(V[:,v])[:tolist]() for v=1:size(V,2) ])
	scale = unit*maximum(p.SIZE(Any[1,2])(submodel))/6
	hpc = p.larModelNumbering(1,1,1)(V,PyObject([VV,EV]),submodel,scale)
	p.VIEW(hpc)
end

# Produce an indexed view of Lar model (W,EW,FW)
function viewLarIndices(W::Array{Float64,2},EW::Array{Int64,2},FW::Array{Int64,2},unit=1.0)
	EV = [EW[:,k] for k=1:size(FW,2)]
	viewLarIndices(W,EV,FW)
end

# Produce an indexed view of Lar model (W,EW,FW)
function viewLarIndices(W::Array{Float64,2},EW::Array{Int64,2},,FW::Array{Array{Int64,1},unit=1.0)
	EW = cols2any(EW)
	function convertData(W,EW)
		# Shifted vertices
		Z = zeros(length(W)+size(W,1))
		for k=1:length(W)
			Z[k + size(W,1)] = W[k]
		end
		Z = reshape(Z,size(W,1),size(W,2)+1)
		# shifted edges
		EZ = Any[[1,1]]
		for k=1:length(EW)
			push!(EZ, EW[k]+1)
		end
		Z,EZ
	end
	V,EV = convertData(W,EW)
	submodel = lar2hpc(V,hcat(EV...))
	VV = PyObject(Any[Any[k] for k in  0:size(V,2)-1])
	EZ = map(Array{Int32},EV-1)
	EV = PyObject(Any[PyObject(EZ[e])[:tolist]() for e=1:length(EZ) ])
	FV = PyObject(Any[PyObject(FW[f])[:tolist]() for f=1:length(FW) ])
	V = PyObject(Any[PyObject(V[:,v])[:tolist]() for v=1:size(V,2) ])
	Scale = unit*maximum(p.SIZE(Any[1,2])(submodel))/6
	hpc = p.larModelNumbering(1,1,1)(V,PyObject([VV,EV,FV]),submodel,Scale)
	p.VIEW(hpc)
end


# end # module Lar-core
