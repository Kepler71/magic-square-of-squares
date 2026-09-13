P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve((-2174406)*x^1+(5026188)*x^2+(-5026188)*x^4+(2174406)*x^5); J:=Jacobian(H);
pts:=[R : R in Points(H : Bound:=2000) | R[3] ne 0 and R[2] ne 0];
inf:=[R : R in Points(H : Bound:=2) | R[3] eq 0][1];
xs:={R[1]/R[3] : R in pts}; print "XS:", xs;
G:=[J!(R-inf) : R in pts | R[2] gt 0];
M:=HeightPairingMatrix(G); print "HPM:", M; print "DET:", Determinant(M);
