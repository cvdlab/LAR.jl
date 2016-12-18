using LAR

v,(vv,ev,fv,cv) = p.larCuboids((2,2,1),true)
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
Z,FZ = cmerge(models)
viewexploded(Z,FZ)

