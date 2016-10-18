
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
function edgeSlopeOrdering(VE,V)
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
	forwardBackward = edgeSlopeOrdering(VE,V)
	for (v,(edges,nexts,prevs)) in enumerate(forwardBackward)
		for (h,e) in enumerate(edges)
			S[edges[h],v] = (nexts[h],prevs[h])
		end
	end
	return S
end

S = cycleBasis(W,EW)

