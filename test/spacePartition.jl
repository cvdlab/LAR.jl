"""
# include("test/merge.jl")

The `spacePartition` function takes as input a **non-valid** (with the meaning used in solid modeling field --- see~\cite{Requicha:1980:RRS:356827.356833}) `LAR` model of dimension $d-1$, i.e.~a triple `(V,FV,EV)`, and an array `buckets` indexed on faces, and containing the subsets of faces with greatest probability of intersecting each indexing face, respectively. 

The `spacePartition` function returns the **valid** `LAR` boundary model `(W,FW,EW)` of the space partition induced by `FV`.
"""
function spacePartition(V,FV,EV, buckets)

	""" input: face index f; candidate incident faces F[f]; """
	
	for (f,F) in enumerate(buckets)
	
		""" Computation of submanifold map M moving f to z=0 """
		
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