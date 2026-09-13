import sys
from chab_all_lib import poly_coeffs, magma
A,C,d3=map(int,sys.argv[1:4])
f=poly_coeffs(d3,A,C)
print(magma(f"""P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve(P!{f});
pts,flag:=RationalPointsGenus2(H);
print "NPTS:", #pts; print "PROVEN:", flag;
print "TVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];
"""))
