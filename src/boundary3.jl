using LAR

include("test/merge.jl")
include("test/spacePartition.jl")



v,(vv,ev,fv,cv) = p.larCuboids((4,4,2),true)
V = hcat([Array{Float64,1}(v[k,:]) for k=1:size(v,1)]...)
FV = hcat([Array{Int64,1}(fv[k,:]+1) for k=1:size(fv,1)]...)
EV = hcat([Array{Int64,1}(ev[k,:]+1) for k=1:size(ev,1)]...)
model1 = Any[V,FV,EV]

W = hcat([V[:,k] + [.5;0.5;0.] for k=1:size(V,2)]...)
W = rotate([0,0,π/6],W)
FW = copy(FV)
EW = copy(EV)
model2 = Any[W,FW,EW]

models = [model1,model2]
X,FX,EX = cmerge(models)

viewexploded(X,FX)
viewexploded(X,EX)



V,FV,EV = deepcopy((X,FX,EX))
W,FW,EW = spacePartition(V,FV,EV)

viewexploded(W,EW)
viewexploded(W,FW)

#FE = crossRelation(FW,EW)
#boundaryTriangulation(W,FW,EW,FE)



