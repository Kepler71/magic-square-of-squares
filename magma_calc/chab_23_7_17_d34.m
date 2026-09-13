P<x>:=PolynomialRing(Rationals());
T0:=Cputime();
H:=HyperellipticCurve(34*x*(x^2-1)*(169*x-409)*(409*x-169));
H2,phi:=ReducedMinimalWeierstrassModel(H);
print "MODEL:", H2;
J:=Jacobian(H2);
Q1:=J![x^2+4*x, 240*x];
Q2:=J![x^2-9/2*x-34, 384/5*x-3264/5];
Q3:=J![x^2+25/2*x+34, 640/3*x+5440/3];
print "ORDERS:", Order(Q1), Order(Q2), Order(Q3), "t:", Cputime(T0);
inf:=[Q : Q in [Q1,Q2,Q3] | Order(Q) eq 0];
if #inf gt 0 then
  pts2:=Chabauty(inf[1]);
  pts:=[Inverse(phi)(R) : R in pts2];
  print "CHABAUTY_NPTS:", #pts; print "PTS:", pts; print "t:", Cputime(T0);
end if;
