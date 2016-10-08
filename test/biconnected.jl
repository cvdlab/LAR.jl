#=	Biconnected components of a 1-complex.
An implementation of the Hopcroft-Tarjan algorithm~\cite{Hopcroft:1973:AEA:362248.362272} for computation of the biconnected components of a graph is given here =#


# Adjacency lists of 1-complex vertices 
function vertices2vertices(model)
    V,EV = model
    csrEV = larcc.csrCreate(EV)
    csrVE = larcc.csrTranspose(csrEV)
    csrVV = larcc.matrixProduct(csrVE,csrEV)    
    cooVV = csrVV.tocoo()
    data,rows,cols = AA(list)([cooVV.data, cooVV.row, cooVV.col])
    triples = zip(data,rows,cols)
    VV = [[] for k in range(len(V))]
    for datum,row,col in triples
        if row != col: VV[col] += [row]
    end
    AA(sorted)(VV)
end

# Main procedure for biconnected components 
function biconnectedComponent(model)
    W,_ = model
    V = range(len(W))
    count = 0
    stack,out = [],[]
    visited = [None for v in V]
    parent = [None for v in V]
    d = [None for v in V]
    low = [None for v in V]
    for u in V: visited[u] = False
    for u in V: parent[u] = []
    VV = vertices2vertices(model)
    for u in V
        if not visited[u]
            DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
        end
    end
    W,[component for component in out if len(component) > 1]
end

# Hopcroft-Tarjan algorithm 
function DFV_visit( VV,out,count,visited,parent,d,low,stack,u )
    visited[u] = True
    count += 1
    d[u] = count
    low[u] = d[u]
    for v in VV[u]
        if not visited[v]
            stack += [(u,v)]
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]
                out += [outputComp(stack,u,v)]
            end
            low[u] = min( low[u], low[v] )
        else
            if not (parent[u]==v) and (d[v] < d[u])
                stack += [(u,v)]
                low[u] = min( low[u], d[v] )
            end
        end
    end
end

# Output of biconnected components 
function outputComp(stack,u,v)
    out = []
    while True
        e = stack.pop()
        out += [list(e)]
        if e == (u,v) 
        	break
        end
    end
    list(set(AA(tuple)(AA(sorted)(out))))
end

