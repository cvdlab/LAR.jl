using LAR

# test coarse surface approximation (r=1, R=3)
W,lar = toroidalsurface(1,3,10,20);
larview(W,lar.FV)
larview(W,lar.EV)
larview(W,lar.VV)
viewexploded(W,lar.FV)
viewexploded(W,lar.EV)
viewexploded(W,lar.VV)



# test non-manifold surface  (r=3, R=3)
V,lar = toroidalsurface(1,1,40,80);
larview(V,lar.FV)
viewexploded(V,lar.FV)


# test coarsest surface approximation (r=1, R=3)
W,lar = toroidalsurface(1,1,3,3);
larview(W,lar.FV)
larview(W,lar.EV)
viewexploded(W,lar.FV)
viewexploded(W,lar.EV)


