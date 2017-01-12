
# first example in biconnected.jl
pol = W,hcat([EW[:,e] if e>0 else reverse(EW[:,-e]) for e in chains[1]]...)

cycle = []
for edge in chains[1]
	if edge>0
		push!(cycle, EW[:,edge])
	elseif edge<0
		push!( cycle, Array{Int64,1}(reverse(EW[:,-edge])) )
	end
end
EV = hcat(cycle...)

view(V,EV)
pointInPolygonClassification(V,EV)([200,300])
