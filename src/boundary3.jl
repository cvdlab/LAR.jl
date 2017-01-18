using LAR

include("test/merge.jl")
include("test/spacePartition.jl")


v,(vv,ev,fv,cv) = p.larCuboids((1,1,1),true)
V = hcat([Array{Float64,1}(v[k,:]) for k=1:size(v,1)]...)
FV = hcat([Array{Int64,1}(fv[k,:]+1) for k=1:size(fv,1)]...)
EV = hcat([Array{Int64,1}(ev[k,:]+1) for k=1:size(ev,1)]...)
model1 = Any[V,FV,EV]

W = hcat([V[:,k] + [.5;0.5;0.] for k=1:size(V,2)]...)
#W = rotate([0,0,π/6],W)
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


function chainCoords(EW,rows)
	edges = hcat([EW[:,e] for e in rows]...)
	for k=1:length(edges)
		v1,v2 = find(edges .== edges[k])
		if (v1 + v2) % 2 == 0
			# swap column elements, for column that contains v2
			col = v2 ÷ 2
			if v2 % 2 != 0
				col +=1
			end 
			edges[2,col],edges[1,col] = edges[1,col],edges[2,col]
		end
	end
	out = Int[]
	for k=1:size(edges,2)
		if edges[1,k] < edges[2,k]
			push!(out,1)
		else
			push!(out,-1)
		end
	end
	out
end

I,J,V = Int[],Int[],Int[]
FE = crossRelation(FW,EW)
for (f,rows) in enumerate(FE)
	@show rows
	append!(I,rows)
	@show [f for k=1:length(rows)]
	append!(J,[f for k=1:length(rows)])
	@show chainCoords(EW,rows)
	append!(V,chainCoords(EW,rows))
end

full(sparse(I,J,V))



#boundaryTriangulation(W,FW,EW,FE)

viewLarIndices(W,EW,FW,2)

