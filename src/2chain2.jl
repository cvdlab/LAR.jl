using LAR
include("src/inters.jl")

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

function theSign(EV,edge,v)
	w = pop!(setdiff(Set(EV[edge]), Set(v)))
	v < w ? 1 : -1
end

# Return the sparse matrix codifying (next,prev) edges for each pair [v,e]
function cycleBasis(V,EV)
	I = Int64[]; 
	J = Int64[]; 
	Vnext, Vprev = Int64[], Int64[] 
	for (e,(i,j)) in enumerate(EV)
		push!(I,i); push!(I,j); 
		push!(J,e); push!(J,e); 
		push!(Vnext,0); push!(Vnext,0); 
		push!(Vprev,0); push!(Vprev,0); 
	end
	Snext = sparse(J,I,Vnext)
	Sprev = sparse(J,I,Vprev)
	rows = rowvals(Snext) 
	VE = Array{Array{Int64,1},1}()
	for k=1:size(Snext,2)
		# extract column
		col = [rows[h] for h in nzrange(Snext,k)]
		push!(VE, col)
	end
	forwardBackward = edgeSlopeOrdering(VE,V,EV)
	for (v,(edges,nexts,prevs)) in enumerate(forwardBackward)
		for (h,e) in enumerate(edges)
			esign = theSign(EV,edges[h],v)
			Snext[edges[h],v] = nexts[h]*esign
			Sprev[edges[h],v] = prevs[h]*esign
		end
	end
	return Snext,Sprev
end

# Return the sparse matrix of the signed ∂–1 operator
function boundary_1(EV)
	we = cellComplex(EV)';
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


datafile = readcsv("test/svg/test1.lines");
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2));
len = length(datafile);
EV = collect(reshape(1:(len÷2), 2,(len÷4)));
lineArray = (V,EV);
W,EW = lines2lar(lineArray);
lar = Lar();
lar.EV = [EW[:,k] for k=1:size(EW,2)];
W,newcells,oldindex = larvalidate(W,lar,10^2);
EW = newcells.EV;
viewLarIndices(W,EW,0.75)


using Iterators

function colProd(op,chain1)
	chain0 = op * chain1
	chain0_cells = []
	I,J,V = findnz(chain0)
	@itr for (i,v) in zip(I,V)
    	push!( chain0_cells, i*v )
    end		
	chain0_cells
end

function rowProd(chain0,op)
	chain1 = chain0' * op
	chain1_cells = []
	I,J,V = findnz(chain1)
	@itr for (j,v) in zip(J,V)
    	push!( chain1_cells, j*v )
    end		
	chain1_cells
end

function extract_2cell(seed, Snext,Sprev, useCounts, B_1, EW)
	chain1 = map(Int64, spzeros(size(B_1,2),1));
	if useCounts[seed] == 0
		chain1[seed,1] = +1
	elseif useCounts[seed] == 1
		chain1[seed,1] = -1
	end
	chain1_cells = [chain1[seed,1] * seed]
	chain0_cells = colProd(B_1, chain1)
	useCounts[seed] += 1
	
	while chain0_cells != []
		chain1_cells = rowProd(chain0, B_1)
	
		next1_cells = [ S[abs(cell1),abs(cell0)] for cell0 in chain0_cells, 
			cell1 in chain1_cells if Set(abs(cell0)) < Set(EW[cell1])]
	
		direction_test = all([ useCounts[cell]<2 for cell in next1_cells ])
		if direction_test==true
			for cell in Set(next1_cells)
				useCounts[cell] += 1
				chain1[cell,1] = 1
				push!(chain1_cells, cell)
			end
		else 
			return []
		end
		chain0 = B_1 * chain1
		cells0,val = findnz(chain0[:,1])
		chain0_cells = [cells0[k] for (k,v) in enumerate(val) if abs(val[k])==1]
	end
	chain1_cells
end
			
Snext,Sprev = cycleBasis(W,EW)
B_1 = boundary_1(EW)
S = Snext
useCounts = [0 for k in 1:length(EW)]
extract_2cell(1, S, useCounts, B_1, EW)
			
