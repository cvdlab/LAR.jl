#include("merge.jl")
#include("boundary.jl")

"""
The `spacePartition` function takes as input a **non-valid** (with the meaning used in solid modeling field --- see~\cite{Requicha:1980:RRS:356827.356833}) `LAR` model of dimension $d-1$, i.e.~a triple `(V,FV,EV)`, 

The `spacePartition` function returns the **valid** `LAR` boundary model `(W,FW,EW)` of the space partition induced by `FV`.

First an array `buckets` indexed on faces is computed. It contains the subsets of faces with greatest probability of intersecting each indexing face, respectively. 
"""


function crossRelation(FV::Array{Int64,2},EV::Array{Int64,2})
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = crossRelation(FW,EV)
end

function crossRelation(FV::Array{Array{Int64,1},1},EV::Array{Int64,2})
	sparseFV = cellComplex(FV)
	sparseEV = cellComplex(EV)
	sparseEF = (sparseFV' * sparseEV)'
	I,J,V = findnz(sparseEF)
	FE = [Int64[] for k=1:length(FV)]
	for (i,j,v) in zip(I,J,V) 
		if v==2
			push!(FE[j],i)
		end
	end
	return FE
end



function subModel(V::Array{Float64,2}, FV::Array{Int64,2}, EV::Array{Int64,2}, 
				bucket::Array{Int64,1},f::Int64,FE::Array{Array{Int64,1},1})
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = subModel(V,FW,EV,bucket,f,FE)
end

function subModel(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
					EV::Array{Int64,2},bucket::Array{Int64,1},f::Int64, FE::Array{Array{Int64,1},1})
	FW = [FV[f] for f in bucket]
	EW = hcat(collect(Set([EV[:,e] for f in bucket for e in FE[f]]))...)

	oldVerts = Set(Int64[])
	for k=1:length(bucket)
		union!(oldVerts, Set(FW[k]))
	end	
	olds = collect(oldVerts)
	vdict = Dict(zip(olds, 1:length(olds)))
	fz = [[vdict[v] for v in FW[k]] for k=1:length(FW)]
	ez = hcat([[vdict[v] for v in EW[:,k]] for k=1:size(EW,2)]...)
	z = hcat([V[:,v] for v in olds]...)
	pivot = find(bucket .== f)[1]
	return z,fz,ez,pivot
end


function visualize(V,FV,EV,f)
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = visualize(V,FW,EV,f)
end

function visualize(Z,FZ::Array{Array{Int64,1},1},EZ,f)
	bucket = lar2hpc(Z,FZ)
	bucketEdges = lar2hpc(Z,EZ)
	params = PyObject(pyeval("list([1.,0.,0.,0.1,  0.,1.,0.,0.1,  
				0.,0.,1.,0.1, 0.,0.,0.,0.1, 100.])"))
	glass = p.MATERIAL(params)
	pivot = p.COLOR(p.RED)(p.JOIN(lar2hpc(Z,[FZ[f]])))
	p.VIEW(p.STRUCT([glass(bucket),bucketEdges,pivot]))
end



function spacePartition(V::Array{Float64,2}, FV::Array{Int64,2}, EV::Array{Int64,2},
						debug=false)
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = spacePartition(V,FW,EV,debug)
end


function spacePartition(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
						EV::Array{Int64,2},debug=false)
	""" input: face index f; candidate incident faces F[f]; """
	boxes = lar2boxes(V,FV)
	buckets = boxBucketing(boxes)
	FE = crossRelation(FV,EV)
	
	for (f,F) in enumerate(buckets)
		@show (f,F)
		""" F[f] submodel extraction from (V,FV,EV) """
		Z,FZ,EZ,pivot = subModel(V,FV,EV,F,f,FE)
		if debug visualize(Z,FZ,EZ,pivot) end
		
		""" Computation of submanifold map M moving f to z=0 """
		M = submanifoldMapping(V,FV,f)
		
		""" Transformation of F(f) by M, giving (W,EW,FW) := M(F(f)) """
		Z1 = vcat(Z,ones((1,size(Z,2))))
		Y = M * Z1
		Z = Y[ 1:3, : ]
		if debug visualize(Z,FZ,EZ,pivot) end
		
	end
end

		""" filtering of EW edges traversing z=0, """
		""" giving EZ edges and incident faces FZEZ """
		
		
		
		""" for each face in FZEZ, computation of the aligned set of points p(z=0) """
		
		""" Remove external vertices """
		
		""" Apply the inverse submanifold transform """
		
		""" Accumulate the submodel parts """
		
	end
	
	""" return the **valid** `LAR` 2-skeleton model `(W,FW,EW)` 
	
end

V,FV,EV = deepcopy((X,FX,EX))

chains = boundary(V,EV)
operator = boundaryOp(EW,chains)
#println(full(operator))
FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]















