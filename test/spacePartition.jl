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


function spacePartition(V::Array{Float64,2}, Array{Int64,2}, EV::Array{Int64,2})
	FW = [FV[:,k] for k=1:size(FV,2)]
	out = spacePartition(V,FW,EV)
end

function spacePartition(V::Array{Float64,2}, FV::Array{Array{Int64,1},1}, 
						EV::Array{Int64,2})

	""" input: face index f; candidate incident faces F[f]; """
	boxes = lar2boxes(V,FV)
	buckets = boxBucketing(boxes)
	FE = crossRelation(FV,EV)
	
	for (f,F) in enumerate(buckets)
		@show (f,F)
		Z,EZ,FZ = boxes3lar(lar2boxes(X,[FX[:,f] for f in buckets[cell]]))
		""" Computation of submanifold map M moving f to z=0 """
		M = submanifoldMapping(V,FV,cell)
		
		
		
		ZZ = vcat(Z,ones((1,size(Z,2))))
		Y = M*ZZ
		bucket = lar2hpc(Y[1:3,1:size(Z,2)],FZ)
		bucketEdges = lar2hpc(Y[1:3,1:size(Z,2)],EZ)
		
		""" Transformation of S(f) by M, giving S = (W,EW) := M(S(f)) """
		
		""" filtering of EW edges traversing z=0, """
		""" giving EZ edges and incident faces FZEZ """
		
		""" for each face in FZEZ, computation of the aligned set of points p(z=0) """
		
		""" Remove external vertices """
		
		""" Apply the inverse submanifold transform """
		
		""" Accumulate the submodel parts """
		
	end
	
	""" return the **valid** `LAR` 2-skeleton model `(W,FW,EW)` 
	
end

V,FV,EV = X,FX,EX

chains = boundary(V,EV)
operator = boundaryOp(EW,chains)
#println(full(operator))
FW = [sort(collect(Set(vcat([EW[:,abs(e)] for e in face]...)))) for face in chains]















