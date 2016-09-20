V = Verts([[1,2] [3,4] [5,6]])
m,n = size(V)
mymodel = Lar()
fieldnames(Lar)
fieldnames(Verts)

V = Float64[[1,2] [3,4] [5,6]]
vs = Verts(V)

mymodel = Lar()
fieldnames(Lar)
fieldnames(Verts)


V = readcsv("/Users/paoluzzi/Documents/dev/mtj/testfile.V",Float64)
EV = readcsv("/Users/paoluzzi/Documents/dev/mtj/testfile.EV",Int32)
EV1 = readcsv("/Users/paoluzzi/Documents/dev/mtj/testfile1.EV")

lar = JSON.parsefile("/Users/paoluzzi/Documents/dev/mtj/larfile.json"; 
		dicttype=Dict, use_mmap=true)

lar = JSON.parse("""
	{"V":[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]],
	"EV":[[1,2],[2,3],[3,4],[4,1]],
	"FV":[[4,1,3,2]]}
""")

larstring = """
	{"V" : [[5,0],[7,1],[9,0],[13,2],[15,4],[17,8],[14,9],[13,10],
	[11,11],[9,10],[5,9],[7,9],[3,8],[0,6],[2,3],[2,1],[8,3],
	[10,2],[13,4],[14,6],[13,7],[12,10],[11,9],[9,7],[7,7],[4,7],
	[2,6],[3,5],[4,2],[6,3],[11,4],[12,6],[12,7],[10,5],[8,5],
	[7,6],[5,5]],

	"FV" : [[0,1,16,28,29],[0,15,28],[1,2,17],[1,16,17,33],[2,3,17],
	[3,4,18,19],[3,17,18,30],[4,5,19],[5,6,19],[6,7,20,21,22,32],
	[6,19,20],[7,8,21],[8,9,21,22],[9,11,23,24],[9,22,23],
	[10,11,24,25],[10,12,25],[12,13,25,26],[13,14,27],[13,26,27],
	[14,15,28],[14,27,28,29,36],[16,29,34],[16,33,34],[17,30,33],
	[18,19,31],[18,30,31],[19,20,31,32],[22,23,32,33],[23,24,34,35],
	[23,33,34],[24,25,27,36],[24,35,36],[25,26,27],[29,34,35],
	[29,35,36],[30,31,32,33]]}
""";

lar = json2larmodel(larstring)
fieldnames(lar)
typeof(lar.Lar)
lar.Lar.VV, lar.Lar.EV, lar.Lar.CV
lar.Lar.FV
map(Array{Int,1}, lar.Lar.FV)


EV1 = readcsv("/Users/paoluzzi/Documents/dev/mtj/testfile1.EV")
EV1 = listOfList(EV1)

model1 = json2larmodel("""
	{"V":[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]],
	"EV":[[1,2],[2,3],[3,4],[4,1]],
	"FV":[[4,1,3,2]]}
""");

mod2 = json2larmodel(larstring);



M_1 = cellComplex(model1.Lar.EV);
M_2 = cellComplex(model1.Lar.FV);
print(full(M_1))
print(full(M_2))


cube = p.COLOR(p.RED)(p.CUBE(0.5))

p.VIEW(cube)

V = Float64[[1,2] [3,4] [5,6]]
a = PyObject(V)

verts = PyObject(a[:tolist]())


println("V =\n",V)
println("EV =\n",EV)
println("EV1 =\n",EV1)

methods(view)

larmodel = p.larCuboids((2,2,2),true);
verts,cells,chainbases = importmodel(larmodel);
view(verts,cells.EV);
view(verts,cells.FV);
view(verts,cells.CV);
view(verts,cells.VV);



v,cv = p.larCuboids((1,2,1))

v,cv = p.larCuboids((20,20,20));
viewexploded(v,cv)

larmodel = p.larCuboids((2,2,2),full=true);
(v,cc)=larmodel;
(vv,ev,fv,cv) = cc;
cc = p.CAT([PyObject(vv),PyObject(ev),PyObject(fv),PyObject(cv)]);
scaleargs=2;
viewexploded(v,cc,scaleargs);




mod1 = json2larmodel("""
	{"V":[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]],
	"EV":[[1,2],[2,3],[3,4],[4,1]],
	"FV":[[4,1,3,2]]}
""")


view(mod1.Verts',mod1.Lar.FV)
view(mod2.Verts',rebase(mod2.Lar.FV))
view(embed(mod2.Verts)',rebase(mod2.Lar.FV))
viewexploded(embed(mod2.Verts)',(mod2.Lar.FV))



larmodel = p.larCuboids((1,1));
(v,fv)=larmodel;
view(translate(v,[0.5 0.5]),rebase(fv))

larmodel = p.larCuboids((2,2,2));
(v,cv)=larmodel;
viewexploded(translate(v,[-2 0 0]),cv,1)


A = reshape(1:30,3,10);
s = [1,1/2,1/3];
scale( A, s );

larmodel = p.larCuboids((2,2,2));
(v,cv)=larmodel;
viewexploded(scale(v,[-1 -1 -1]),cv,1)



larmodel = p.larCuboids((1,1,1));
(v,cv)=larmodel;
viewexploded(rotate([1.,1,1],v')',cv);



