using LAR

include("test/merge.jl")
include("test/spacePartition.jl")

using PyCall
@pyimport triangle


v,(vv,ev,fv,cv) = p.larCuboids((2,2,2),true)
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

function boundary_2(EW,FW)
	I,J,V = Int[],Int[],Int[]
	FE = crossRelation(FW,EW)
	for (f,rows) in enumerate(FE)
		append!(I,rows)
		append!(J,[f for k=1:length(rows)])
		append!(V,chainCoords(EW,rows))
	end
	return sparse(I,J,V)
end

full(boundary_2(EW,FW))

viewLarIndices(W,EW,FW,2)

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


function larTriangulate(W::Array{Float64,2},FW::Array{Int64,2},EW::Array{Int,2})
	FV = [face for face in FW]
	return larTriangulate(W,FV,EW)
end

function larTriangulate(W::Array{Float64,2},FW::Array{Array{Int64,1},1},EW::Array{Int,2})
	FE = crossRelation(FW,EW)
	TW = Array{Array{Int64,2},1}()
	for (f,face) in enumerate(FW)
		# submanifold mapping and 2D projection
		vs = hcat([W[:,v] for v in face]...)
		fv = [[v for v=1:length(face)]]
		ev = Array{Array{Int64,1},1}([])
		for e in FE[f] 
			edge = Int[]
			push!( edge, find(face .== EW[1,e])[1] ) 
			push!( edge, find(face .== EW[2,e])[1] ) 
			push!(ev,edge)
		end
		ev = hcat(ev ...)
		M = submanifoldMapping(vs,fv,1)
		U = vcat(vs,ones((1,size(vs,2))))
		Y = M * U
		Z = Y[ 1:2, : ]
		# prepare data for triangle.triangulation
		polygon = Dict(
			#"holes" => Array{Float64,2}(),
			"segments" => Array{Int32,2}(ev-1)',
			#"vertex_attributes" => hcat([[1.0] for k=1:size(Z,2)]...)',
			"vertices" => Array{Float64,2}(Z'),
		)
		tria = triangle.triangulate(PyObject(polygon), PyObject("p"))
		triangles = [face[v+1] for v in tria["triangles"]]
		push!(TW,triangles)
	end
	return TW
end

function boundary_3(W,EW,FW)
	BND_2 = boundary_2(EW,FW)
	m = size(BND_2,1)
	n = size(BND_2,2)
	EF = crossRelation(EW,FW)

	#forwardBackward = edgeSlopeOrdering(VE,V,EV)
	
end

TW = larTriangulate(W,FW,EW)
viewexploded(W,vcat(TW...)')
