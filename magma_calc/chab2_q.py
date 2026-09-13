import sys
from chab_all_lib import poly_coeffs, magma
A,C,d3=map(int,sys.argv[1:4])
f=poly_coeffs(d3,A,C)
polyF="+".join(f"({c})*x^{i}" for i,c in enumerate(f) if c)
code=f"""P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve({polyF});
J:=Jacobian(H);
pts:=[R : R in Points(H : Bound:=2000) | R[3] ne 0 and R[2] ne 0];
inf:=[R : R in Points(H : Bound:=2) | R[3] eq 0][1];
ptJ:=J!(pts[1]-inf);
print "ORDER:", Order(ptJ);
if Order(ptJ) eq 0 then S:=Chabauty(ptJ); print "NPTS:", #S; print "TVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; end if;
"""
print(magma(code))
