
# vectorized parametric toroidal surface
function toroidalmap(V,r,R,m=40,n=80,angle1=2pi,angle2=2pi)
	V = map(Float64,V)
	V = (scale(V, [angle1/m angle2/n]))
	x =  (R .+ r.*cos(V[:,1])) .* cos(V[:,2])
	y =  (R .+ r.*cos(V[:,1])) .* sin(V[:,2])
	z = -r .* sin(V[:,1])
	[x y z]'
end


function toroidalsurface(r,R,m=10,n=20,back=1,angle1=2pi,angle2=2pi)
	larmodel = p.larCuboids((m,n),true)
	verts,cells,chainbases = importmodel(larmodel)
	V = toroidalmap(verts',r,R,m,n)
	if back==1 
		return V,cells
	elseif 1 < back < 4	
		W,newcells,oldindex = refactoring_incidence(V,cells)
		if back==2
			return W,newcells
		else # back = 3
			Z = map(Float64,hcat([V[:,k] for k in oldindex]...))
			return Z,newcells
		end
	else print("Error: wrong return in `toroidalsurface` ")
	end
end

function refactoring_incidence(V,lar)
	W = zeros!(map(round, V*10^5)/10^5)
	vdict = Dict([ ("$(W[:,k])",k) for k=1:size(W,2) ])
	verts = Dict(zip( keys(vdict), 1:length(vdict) ))
	newindex = convertindex(W,verts)
	oldindex = invertindex(newindex)
	vs = [eval(parse(key))' for key in keys(vdict)]
	W = vcat(vs...)'
	newcells = Lar()
	FW = relink(lar.FV,newindex)
	newcells.FV = [ FW[:,k][:] for k=1:size(FW,2) ]
	W,newcells,oldindex
end

function relink(basis,newindex)
	m0,n0 = length(basis),length(basis[1])
	FV = reshape([basis...],n0,m0)
	FW = round(Int64, zeros(size(FV)))
	for k in 1:size(FV,2)
		for h in 1:size(FV,1)
			FW[h,k] = newindex[FV[h,k]]
		end
	end
	FW
end



V,lar = toroidalsurface(1,3,10,20)
view(V,lar.FV)
viewexploded(V,lar.FV)



W,lar = toroidalsurface(1,3,10,20,3)
view(W,lar.FV)
viewexploded(W,lar.FV)


W,lar = toroidalsurface(3,3,10,20,3)
view(W,lar.FV)
viewexploded(W,lar.FV)


