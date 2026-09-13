P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1-3*z)*(1-86^2*z^2)*(1-89^2*z^2));
pts,flag:=RationalPointsGenus2(C);
print "PROVEN:", flag; print "NPTS:", #pts;
print "ZVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];
