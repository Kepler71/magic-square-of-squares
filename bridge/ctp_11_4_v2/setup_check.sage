from sage.all import *
m,n=11,4
s=QQ(137)/2
b=s*m*m*n*n
print("b =",b)
roots=[-b,-s*m**4,-s*n**4]
print("roots(E orig) =",roots)
R=PolynomialRing(QQ,'x');x=R.gen()
f=prod(x-e for e in roots)
print("f =",f)
E=EllipticCurve([0,f[2],0,f[1],f[0]])
print("E =",E)
M=E.minimal_model()
print("M =",M)
print("M ainvs",M.ainvs())
iso=E.isomorphism_to(M)
print("iso",iso)
mr=[iso(E([e,0]))[0] for e in roots]
print("minimal roots",mr)
print("sum",sum(mr))
print("c4,c6 of M:",M.c4(),M.c6())
I=M.c4();J=2*M.c6()
print("I=",I);print("J=",J)
print("disc16/27:",16*(4*I**3-J**2)/27)
print("cond",M.conductor(),factor(M.conductor()))
print("disc M",M.discriminant(),factor(M.discriminant()))
