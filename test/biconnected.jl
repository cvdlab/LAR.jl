#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#

using LAR
include("src/inters.jl")


# Adjacency lists of 1-complex vertices 
function vertices2vertices(V,EV)
	lar = Lar()
	lar.EV = [EV[:,k] for k=1:size(EV,2)]
	ev = cellComplex(lar.EV)
	vv = ev'*ev
	triples = hcat(findnz(vv)...)'
	VV = [Int64[] for k=1:size(vv,1)]
	for k=1:size(triples,2)
		row,col,datum = triples[:,k]
		if row != col
			push!(VV[col],row)
		end
	end
	return VV
end

# Main procedure for biconnected components 
function biconnectedComponent(W,EV)
    V = 1:length(W)
    count = 0
    stack,out = [],[]
    visited = [false for v in V]
    parent = [-1 for v in V]
    d = [-1 for v in V]
    low = [-1 for v in V]
    VV = vertices2vertices(W,EV)
    for u in V
        if !visited[u]
            DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
        end
    end
    out = [component for component in out if length(component) > 1]
    EV = [[u,v] for (u,v) in vcat(out...)]
    hcat(EV...)
end



# Hopcroft-Tarjan algorithm 
function DFV_visit( VV,out,count,visited,parent,d,low,stack,u )
    visited[u] = true
    count += 1
    d[u] = count
    low[u] = d[u]
    if u > length(VV) return() end
    for v in VV[u] 
        if !visited[v]
            push!(stack, (u,v))   # possible bug: []
            #println("stack = $stack")
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]
                push!(out,outputComp(stack,u,v))
                #println("\tout = $out")
                #println("\tv = $v")
                if v == length(VV) return() end
            end
            low[u] = min( low[u], low[v] )
        else
            if ! (parent[u]==v) & (d[v] < d[u])
                push!(stack, (u,v))   # possible bug: []
                low[u] = min( low[u], d[v] )
            end
        end
    end
end

# Output of biconnected components 
function outputComp(stack,u,v)
    out = []
    while true
        elem = pop!(stack)
        push!(out, elem)
        if elem == (u,v) 
        	break
        end
    end
    out
end


datafile = readcsv("test/svg/test1.lines")
V = reshape(datafile',(size(datafile',1)÷2,size(datafile',2)*2))
len = length(datafile)
EV = collect(reshape(1:(len÷2), 2,(len÷4)))
view(V,EV)
lineArray = (V,EV)
W,EW = lines2lar(lineArray)
lar = Lar()
lar.EV = [EW[:,k] for k=1:size(EW,2)]
W,newcells,oldindex = larvalidate(W,lar)
EW = newcells.EV
viewexploded(W,EW)


lineArray = randomLines(300,.3);
V,EV = lines2lar(lineArray);
viewexploded(V,EV)
VV = vertices2vertices(V,EV)
EW = biconnectedComponent(V,EV)
viewexploded(V,EW)

