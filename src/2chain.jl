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

S = cycleBasis(EW)

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

