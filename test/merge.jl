using LAR
using IntervalTrees
include("src/inters.jl")

v,(vv,ev,fv,cv) = p.larCuboids((10,10,1),true)
V = hcat([Array{Float64,1}(v[k,:]) for k=1:size(v,1)]...)
FV = hcat([Array{Int64,1}(fv[k,:]+1) for k=1:size(fv,1)]...)
model1 = Any[V,FV]

W = hcat([V[:,k] + [.5;.5;.5] for k=1:size(V,2)]...)
FW = copy(FV)
model2 = Any[W,FW]

function cmerge(models)
	V = hcat([models[k][1] for k=1:length(models)]...)
	shifts = [0]
	append!(shifts, [size(models[h][1],2) for h in 1:length(models)])
	FV = hcat([models[k][2]+shifts[k] for k=1:length(models)]...)
	return V,FV
end

models = [model1,model2]
X,FX = cmerge(models)
viewexploded(X,FX)

boxes = lar2boxes(X,FX)
buckets = boxBucketing(boxes)

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


params = PyObject(pyeval("list([1.,0.,0.,0.1,  0.,1.,0.,0.1,  0.,0.,1.,0.1, 0.,0.,0.,0.1, 100.])"))
glass = p.MATERIAL(params)

cell = 111
Z,EZ,FZ = boxes3lar(lar2boxes(X,[FX[:,f] for f in buckets[cell]]))
pivot = p.COLOR(p.RED)(p.JOIN(lar2hpc(boxes3lar(boxes[:,cell])[1:2]...)))
bucket = lar2hpc(Z,FZ)
p.VIEW(p.STRUCT([glass(bucket),pivot]))
