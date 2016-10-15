function cycleBasis(EV)
	I = Int64[]; 
	J = Int64[]; 
	V = Int64[]; #Tuple{Int64,Int64}[] 
	for (e,(i,j)) in enumerate(EV)
		push!(I,i); 
		push!(I,j); 
		push!(J,e); 
		push!(J,e); 
		push!(V,1); #push!(V,(0,0));   
		push!(V,1); #push!(V,(0,0)); 
	end
	S = sparse(I,J,V)
end

S = cycleBasis(EV)

S = S'
rows = rowvals(S) 
for k=1:size(S,2)
	print(k," ")
	# extract column
	for h in nzrange(S,k)
		print(rows[h]," ")
	end
	println()
end


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

viewLarIndices(W,EW)
