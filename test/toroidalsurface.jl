using LAR

# test coarse surface approximation (r=1, R=3)
W,lar = toroidalsurface(1,3,10,20);
view(W,lar.FV)
view(W,lar.EV)
view(W,lar.VV)
viewexploded(W,lar.FV)
viewexploded(W,lar.EV)
viewexploded(W,lar.VV)



# test non-manifold surface  (r=3, R=3)
V,lar = toroidalsurface(1,1,40,80);
view(V,lar.FV)
viewexploded(V,lar.FV)


# test coarsest surface approximation (r=1, R=3)
W,lar = toroidalsurface(1,1,3,3);
view(W,lar.FV)
view(W,lar.EV)
viewexploded(W,lar.FV)
viewexploded(W,lar.EV)


