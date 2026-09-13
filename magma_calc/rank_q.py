import sys,json
from chab_all_lib import poly_coeffs, magma
A,C,d3=map(int,sys.argv[1:4]); extra=sys.argv[4] if len(sys.argv)>4 else ''
f=poly_coeffs(d3,A,C)
polyF="+".join(f"({c})*x^{i}" for i,c in enumerate(f) if c)
code=f"""P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve({polyF});
J:=Jacobian(H);
print "PTS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in Points(H : Bound:=2000)];
print "RANKBOUND:", RankBound(J);
{extra}
"""
print(magma(code))
