ComplexMerge story :-)
======================

I started to see the true story of our algorithm.

1)	It was made possible by a change of language
---------------------------------------------

We left the standard language of  computer science, based on algorithms and data structures, using ONLY standard mathematical language of algebraic topology: space decompositions, cellular complexes, (co)chain linear spaces, boundary and coboundary as linear operators.

Nothing more.


2)	Tooling consisted in (always) lowering the dimension:

In order to merge two 3-complexes, we start from the single bunch of their 2-cells.

for each 2-cell alpha
1.	compute their discrete neighbourhood 
2.	roto-translate the space, moving alpha (and incident 2-cells) to z=0 (considering local submanifold in yje topology of the 2-skeleton?)
3.	consider the boundaries of such minimal problem (1-cells)
4.	intersect such lines in 2D, creating a 2D-partition of alpha
