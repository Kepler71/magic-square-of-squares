import sys
from chab_all_lib import poly_coeffs, magma
A,C,d3=map(int,sys.argv[1:4])
f=poly_coeffs(d3,A,C)
polyF="+".join(f"({c})*x^{i}" for i,c in enumerate(f) if c)
print(magma(f"""P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve({polyF}); J:=Jacobian(H);
pts:=[R : R in Points(H : Bound:=2000) | R[3] ne 0 and R[2] ne 0];
inf:=[R : R in Points(H : Bound:=2) | R[3] eq 0][1];
xs:={{R[1]/R[3] : R in pts}}; print "XS:", xs;
G:=[J!(R-inf) : R in pts | R[2] gt 0];
M:=HeightPairingMatrix(G); print "HPM:", M; print "DET:", Determinant(M);
"""))
