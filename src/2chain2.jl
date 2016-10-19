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

# Return the sparse matrix codifying (next,prev) edges for each pair [v,e]
function cycleBasis(V,EV)
	I = Int64[]; 
	J = Int64[]; 
	Val = Tuple{Int64,Int64}[] 
	for (e,(i,j)) in enumerate(EV)
		push!(I,i); push!(I,j); 
		push!(J,e); push!(J,e); 
		push!(Val,(0,0)); push!(Val,(0,0)); 
	end
	S = sparse(J,I,Val)
	rows = rowvals(S) 
	VE = Array{Array{Int64,1},1}()
	for k=1:size(S,2)
		# extract column
		col = [rows[h] for h in nzrange(S,k)]
		push!(VE, col)
	end
	forwardBackward = edgeSlopeOrdering(VE,V,EV)
	for (v,(edges,nexts,prevs)) in enumerate(forwardBackward)
		for (h,e) in enumerate(edges)
			S[edges[h],v] = (nexts[h],prevs[h])
		end
	end
	return S
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

S = cycleBasis(W,EW)
B_1 = boundary_1(EW)
