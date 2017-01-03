using LAR
include("src/inters.jl")
include("test/boundary.jl")

function boxes3lar(boxes)
	V = Array{Float64,1}[]
	EV = Array{Int64,1}[]
	FV = Array{Int64,1}[]
	for k=1:size(boxes,2)
		xm,ym,zm,xM,yM,zM = boxes[:,k]
		push!(V,[xm,ym,zm])
		push!(V,[xM,yM,zm])
		push!(V,[xm,yM,zm])
		push!(V,[xM,ym,zm])
		push!(V,[xm,ym,zM])
		push!(V,[xM,yM,zM])
		push!(V,[xm,yM,zM])
		push!(V,[xM,ym,zM])

		push!(EV,[ 8(k-1)+1, 8(k-1)+3 ])
		push!(EV,[ 8(k-1)+1, 8(k-1)+4 ])
		push!(EV,[ 8(k-1)+2, 8(k-1)+3 ])
		push!(EV,[ 8(k-1)+2, 8(k-1)+4 ])
		push!(EV,[ 8(k-1)+4+1, 8(k-1)+4+3 ])
		push!(EV,[ 8(k-1)+4+1, 8(k-1)+4+4 ])
		push!(EV,[ 8(k-1)+4+2, 8(k-1)+4+3 ])
		push!(EV,[ 8(k-1)+4+2, 8(k-1)+4+4 ])
		push!(EV,[ 8(k-1)+1, 8(k-1)+4+1 ])
		push!(EV,[ 8(k-1)+2, 8(k-1)+4+2 ])
		push!(EV,[ 8(k-1)+3, 8(k-1)+4+3 ])
		push!(EV,[ 8(k-1)+4, 8(k-1)+4+4 ])
		
		push!(FV,[ 8(k-1)+1, 8(k-1)+2, 8(k-1)+3, 8(k-1)+4])
		push!(FV,[ 8(k-1)+4+1, 8(k-1)+4+2, 8(k-1)+4+3, 8(k-1)+4+4])
		push!(FV,[ 8(k-1)+1, 8(k-1)+3, 8(k-1)+4+1, 8(k-1)+4+3])
		push!(FV,[ 8(k-1)+1, 8(k-1)+4, 8(k-1)+4+1, 8(k-1)+4+4])
		push!(FV,[ 8(k-1)+2, 8(k-1)+3, 8(k-1)+4+2, 8(k-1)+4+3])
		push!(FV,[ 8(k-1)+2, 8(k-1)+4, 8(k-1)+4+2, 8(k-1)+4+4])
	end
	return hcat(V...), hcat(EV...), hcat(FV...)
end

function boxes2lar(boxes)
	V = Array{Float64,1}[]
	EV = Array{Int64,1}[]
	FV = Array{Int64,1}[]
	for k=1:size(boxes,2)
		xm,ym,xM,yM = boxes[:,k]
		push!(V,[xm,ym])
		push!(V,[xM,yM])
		push!(V,[xm,yM])
		push!(V,[xM,ym])
		push!(EV,[ 4(k-1)+1, 4(k-1)+3 ])
		push!(EV,[ 4(k-1)+1, 4(k-1)+4 ])
		push!(EV,[ 4(k-1)+2, 4(k-1)+3 ])
		push!(EV,[ 4(k-1)+2, 4(k-1)+4 ])
		push!(FV,[ 4(k-1)+1, 4(k-1)+2, 4(k-1)+3, 4(k-1)+4])
	end
	return hcat(V...), hcat(EV...), hcat(FV...)
end


function splitting(boxes,nboxes,ncoords,bucket, finalBuckets,splittingStack,batch=14)
	k, theBucket = bucket
	k += 1	
	k = k % ncoords == 0 ? ncoords : k % ncoords
	threshold = ( sum(boxes[k,:])+sum(boxes[k+ncoords,:]) )/(2*nboxes)
	@show threshold
	pre,into,post = [],[],[]
	for h in theBucket
		themin, themax = boxes[k,h], boxes[k+ncoords,h]
		if themin <= threshold <= themax  
			push!(into, h)
		elseif themax < threshold 
			push!(pre, h)
		elseif themin > threshold 
			push!(post, h)
		end
	end
	current = Set(theBucket) 
	@show p.AA(length)([pre,into,post])
	
	appended = 0
	if length(pre) <= batch
		push!(pre, copy(into))
		append!(finalBuckets, pre)
		appended += 1
	end
	if length(post) <= batch
		push!(post, copy(into))
		append!(finalBuckets, post)
		appended += 2
	end
	if appended==1 
		append!(post, copy(into))
		newblock = (k,post)
		if noloop(splittingStack, ncoords, newblock,current)
			if newblock <= batch
				push!(splittingStack, newblock)
			else 
				push!(finalBuckets, newblock[2])
			end
		else
			push!(finalBuckets, newblock[2])
		end
		return
	elseif appended==2
		append!(pre, copy(into))
		newblock = (k,pre)
		if noloop(splittingStack, ncoords, newblock,current)
			if newblock <= batch
				push!(splittingStack, newblock)
			else 
				push!(finalBuckets, newblock[2])
			end
		else
			push!(finalBuckets, newblock[2])
		end
		return
	elseif appended==3
		return
	elseif appended==0 
		append!(pre, copy(into))
		push!(splittingStack, (k,pre))
		append!(post, copy(into))
		push!(splittingStack, (k,post))
	end	
end

function noloop(splittingStack, ncoords, newblock,current)
	j,newterm = 0,Set(newblock[2])
	if newterm == current
		return false
	end
	len = length(splittingStack)
	while j < ncoords
		stored = Set(splittingStack[len-j][2])
		if stored == newterm 
			return false
		end
		j += 1
	end
	return true
end

function boxIndexing(boxes)
	m = size(boxes,1)
	nboxes = size(boxes,2)
	index = [k for k=1:nboxes]
	ncoords = Int(m/2)

    bucket = (2,collect(1:size(boxes,2)))
    splittingStack = Tuple{Int64,Array{Any,1}}[bucket]
    finalBuckets = Any[]
    while splittingStack != []
        bucket = pop!(splittingStack)
		splitting(boxes,nboxes,ncoords,bucket, finalBuckets,splittingStack)
	end
	return finalBuckets
end


function boxBuckets2(boxes)
	boxinstances = [ Set{Int64}([]) for k=1:size(boxes,2) ]
	finalBuckets = boxIndexing(boxes)
	for bucket in finalBuckets
		println(bucket)
		for b in bucket
			println(b)
			boxinstances[b] = union(boxinstances[b], Set{Int64}(bucket))
		end
	end
	return [sort(collect(item)) for item in boxinstances]
end


bb = boxBuckets2(boxes)

finalBuckets = boxIndexing(boxes)

boxgroups = [hcat([boxes[:,b] for b in bucket]...) for bucket in finalBuckets]
colors = [p.RED,p.GREEN,p.BLUE,p.CYAN,p.MAGENTA,p.YELLOW,p.BLACK,p.WHITE,p.GRAY,p.BROWN,p.ORANGE,p.PURPLE]
p.VIEW(p.STRUCT([ p.COLOR(colors[k%12+1])(lar2hpc(boxes3lar(group)[1:2]...)) for (k,group) in enumerate(boxgroups) ]))

pairs = collect(zip(boxBuckets(boxes),boxBuckets2(boxes)))

p.AA(p.AA(length))(pairs)

for pair in pairs println(p.AA(length)(pair)) end














