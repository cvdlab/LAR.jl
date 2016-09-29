using LAR

V = rand(2,10000)
V = [V[:,k] for k=1:size(V,2)]
p.VIEW(p.STRUCT(p.AA(p.MK)(V)))