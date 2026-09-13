P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve(-7*x*(x^2-1)*(457*x-793)*(793*x-457));
S:=Points(H : Bound:=1000);
print "NPTS:", #S;
print "TVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S];
J:=Jacobian(H);
lo,hi:=RankBounds(J);
print "RANK:", lo, hi;
