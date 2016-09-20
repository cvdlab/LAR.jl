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

# the `Basis` type as alias for Sparse binay Matrix by Columns
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
function rebase(faces::Array{Any,1})
	Any[Int[v+1 for v in face] for face in faces]
end

# transform a 0-based array to a 1-based array
function rebase(faces::Array{Any,2})
	faces + 1
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
		verts,cells,chainbases
	end
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Int32,2})
	a,b = PyObject(V), PyObject(EV)
	verts = PyObject(a[:tolist]())
	cells = PyObject(b[:tolist]())
	p.VIEW(p.MKPOL([verts,cells,1]))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Any,2}, EV::Array{Any,1})
	a,b = PyObject(V), PyObject(EV)
	p.VIEW(p.MKPOL([a,b,1]))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Float64,2}, EV::Array{Any,1})
	a,b = PyObject(V), PyObject(EV)
	verts = PyObject(a[:tolist]())
	cells = b
	p.VIEW(p.MKPOL([verts,cells,1]))
end

# visualize an HPC value from a Julia pair (Verts,Cells)
function view(V::Array{Any,2}, EV::Array{Any,2})
	a,b = PyObject(V), PyObject(EV)
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
function viewexploded(v::Array{Float64,2}, fv::Array{Any,1}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	a,b = PyObject(v), PyObject(fv)
	verts = PyObject(a[:tolist]())
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((verts,b))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Any,2}, fv::Array{Any,1}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	a,b = PyObject(v), PyObject(fv)
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((a,b))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Float64,2}, fv::Array{Any,2}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	a,b = PyObject(v), PyObject(fv')
	verts = PyObject(a[:tolist]())
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((verts,Array[b]-1))))
end

# visualise an exploded `larlib` pair from a Julia pair (Verts,Cells)
function viewexploded(v::Array{Any,2}, fv::Array{Any,2}, scaleargs=(1.2,))
	sx, sy, sz = scalingargs(scaleargs)
	cells = PyObject(map(Int16,fv))[:tolist]()
	a,b = PyObject(v), PyObject(cells)
	p.VIEW(p.EXPLODE(sx,sy,sz)(p.MKPOLS((a,b))))
end

# embed vertices array in a higher dimensional space, adding `ndims` zero coords
function embed(verts::Verts; ndims=1)
	vs = copy(verts)
	println(typeof(vs))
	thetype = typeof(vs[1,1])
	added = zeros(thetype, (ndims, size(vs)[2]))
	vs = [vs ; added]
	out = Verts(vs)
end

translate the columns of `V` matrix by sum with `t` vector
function translate( t, V )
	broadcast(+,t,V)
end
function translate( V, t )
	broadcast(+,t,V)
end

scale the columns of `V` matrix by product times `s` vector
function scale( s, V )
	broadcast(*,s,V)
end
function scale( V, s )
	broadcast(*,V,s)
end

rotate the columns of `V` matrix by properly using the `args` parameters
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

# end # module Lar-core
