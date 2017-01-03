using IntervalTrees


function boxBucketing(boxes::Array{Float64,2})
	nboxes = size(boxes,2)
	dim = Int(size(boxes,1)/2)
	trees = Array{Any,1}()
	for d=1:dim
		intervals = sort([ (boxes[d,k], boxes[d+dim,k], k) for k=1:nboxes ]);
		tree1D = IntervalTree{Float64, IntervalValue{Float64, Int64}}();
		for triple in intervals
			push!(tree1D, IntervalValue{Float64, Int64}(triple...));
		end
		push!(trees, tree1D)
	end
	iterator = Array{IntervalTrees.IntersectionIterator,1}()
	out = Any[]
	for d=1:dim-1
		push!(iterator, intersect(trees[d],trees[d+1]))

		x_sets = [ Set(Int64[]) for k=1:nboxes ]
		y_sets = [ Set(Int64[]) for k=1:nboxes ]
		xy_sets = [ Set(Int64[]) for k=1:nboxes ]

		for (item_i,item_j) in iterator[end]
			push!(x_sets[item_i.value], item_j.value)
			push!(y_sets[item_j.value], item_i.value)
		end
		for k=1:nboxes
			xy_sets[k] = intersect(x_sets[k],y_sets[k])
		end
		push!(out,xy_sets)
	end
	if length(out)==1
		buckets = out[1]
	elseif length(out)==2
		buckets = [intersect(cellsxy,cellsyz) for (cellsxy,cellsyz) in zip(out[1],out[2])]
	end
	return [sort(collect(bucket)) for bucket in buckets]
end

buckets = boxBucketing(boxes)
