# Circular ordering of edges around vertices
function edgeSlopeOrdering(VE,V,EV)
	VE_sorted = Array{Array{Int64,1},1}()
	forwardBackward = Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}}[]
	for (v,ve) in enumerate(VE)
		ve_angle = []
		if ve != []
			for edge in ve
				w = pop!(setdiff(Set(EV[edge]), Set(v)))
				(x,y) = V[:,w] - V[:,v]
				angle = atan2(y,x)
				push!(ve_angle, 180*angle/pi)
			end
		end
		pairs = sort(collect(zip(ve_angle,ve)))
		sortedEdges = [pair[2] for pair in pairs]
		push!(VE_sorted, sortedEdges)
		theNexts,thePrevs = nextPrev(sortedEdges)
		push!(forwardBackward, (sortedEdges,theNexts,thePrevs))
		
	end
	forwardBackward
end

function boundary_1(EV)
	we = cellComplex(EV);
	rows = rowvals(we);
	for e in 1:size(we,2)
		i,j = collect(nzrange(we,e))
		h,k = rows[i], rows[j]	
		if h < k
			we[h,e] = -1
		else 
			we[k,e] = -1
		end
	end;
	we
end

# Compute the codomain of functions (next,prev) for the edge indices
function nextPrev(edges::Array{Int64,1})
	theNexts,thePrevs = Array{Int64,1}(), Array{Int64,1}()
	n = length(edges)
	for k in 1:n
		push!(theNexts, k==n ? edges[1] : edges[k+1] )
		push!(thePrevs, k==1 ? edges[n] : edges[k-1] )
	end
	theNexts,thePrevs
end

function rowColIncidence(BND_1,n)
	vs,es = findn(BND_1)
	VE = [Int64[] for k=1:n]
	for k in 1:length(vs)
		push!(VE[vs[k]], es[k])
	end
	return VE
end

function edgeSlopeOrdering(VE,V,EV)
 	VE_sorted = Array{Array{Int64,1},1}()
 	forwardBackward = Tuple{Array{Int64,1},Dict{Int64,Int64},Dict{Int64,Int64}}[]
 	for (v,ve) in enumerate(VE)
 		ve_angle = []
 		if ve != []
 			for edge in ve
 				w = pop!(setdiff(Set(EV[:,edge]), Set(v)))
 				(x,y) = V[:,w] - V[:,v]
 				angle = atan2(y,x)
 				push!(ve_angle, 180*angle/pi)
 			end
 		end
 		pairs = sort(collect(zip(ve_angle,ve)))
 		sortedEdges = [pair[2] for pair in pairs]
 		push!(VE_sorted, sortedEdges)
 		theNexts,thePrevs = nextPrev(sortedEdges)
 		next = Dict(zip(sortedEdges,theNexts))
 		prev = Dict(zip(sortedEdges,thePrevs))
 		push!(forwardBackward, (sortedEdges,next,prev))
 	end
 	forwardBackward
end


function checkOrient(BND_1,m,Cell_0,pivot,next)
	sign = Cell_0[1,pivot]
	C_1 = sparse([next],[1],[1],m,1)
	C_0 = (BND_1 * C_1)'
	chain_0 = findnz(C_0)
	_,cells_0,vals_0 = chain_0
	k = findfirst((x->x==pivot),cells_0)
	if vals_0[k]*sign==sign
		return -sign
	else 
		return sign
	end
end

function cellTracking(BND_1,m,n,forwardBackward, V,EV, thecell,thesign)
	chain = [thecell]
	signs = [thesign]
	signedChain = [thecell*thesign]
	C_1 = sparse(chain,[1],signs,n,1)
	C_0 = (BND_1 * C_1)'
	chain_0 = findnz(C_0)
	_,cells_0,vals_0 = chain_0
	while cells_0 != []
		for k in 1:length(cells_0)
			pivot = cells_0[k]
			Cell_0 = sparse([1],[pivot],[vals_0[k]],1,m)
			C_1 = (Cell_0 * BND_1)
			chain_1 = findnz(C_1)
			_,cells_1,vals_1 = chain_1
			hinge = pop!(intersect(Set(cells_1), Set(chain)))
			if Cell_0[1,pivot] == -1
				next = forwardBackward[pivot][3][hinge]
			else
				next = forwardBackward[pivot][2][hinge] end
			thesign = checkOrient(BND_1,n,Cell_0,pivot,next)
			signedCell = Int64(thesign*next)
			if !(signedCell in signedChain)
				signedChain = push!(signedChain,signedCell) end
		end
		chain = [abs(cell) for cell in signedChain]
		signs = [sign(cell) for cell in signedChain]
		C_1 = sparse(chain,[1 for cell in chain],signs,n,1)
		C_0 = (BND_1 * C_1)'
		chain_0 = findnz(C_0)
		_,cells_0,vals_0 = chain_0
	end
	
	return signedChain
end


function boundary(V::Array{Float64,2},EV::Array{Int64,2})
	BND_1 = boundary_1(EV)
	m = size(BND_1,1)
	n = size(BND_1,2)
	VE = rowColIncidence(BND_1,n)
	forwardBackward = edgeSlopeOrdering(VE,V,EV)

	out = Array{Int64,1}[]
	store1 = [100 for k=1:size(EV,2)]
	store2 = [100 for k=1:size(EV,2)]
	cell = 1
	thesign = 1
	while cell != 0
		#@show cell
		
		signedChain = cellTracking(BND_1,m,n,forwardBackward,V,EV, cell,thesign)
		push!(out,signedChain)
		# store the cells in signedChain
		for c in signedChain		
			if store1[abs(c)] == 100
				store1[abs(c)] = sign(c)
			elseif store2[abs(c)] == 100
				store2[abs(c)] = sign(c)
			end
		end
		# get the next seed
		pairs = collect(zip(store1,store2))
		cell = findfirst(c->99<=sum(c)<=101, pairs)
		if cell != 0
			thesign = -(sum(pairs[cell])-100)
		else
			cell = findfirst(c->sum(c)==200, pairs)
			if cell==0 break end
		end
	end
	pairs = collect(zip(store1,store2))
	out
end


function boundaryOp(EV::Array{Int64,2}, chains::Array{Array{Int64,1},1})
	I,J,V = Int64[],Int64[],Int64[]
	for (k,chain) in enumerate(chains)
		append!(I, abs(chain))
		append!(J, [k for row in chain])
		append!(V, sign(chain))
	end
	return sparse(I,J,V,size(EV,2),length(chains))
end















