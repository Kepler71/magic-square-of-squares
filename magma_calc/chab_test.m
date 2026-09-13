P<t>:=PolynomialRing(Rationals());
T0:=Cputime();
H:=HyperellipticCurve(34*t*(t^2-1)*(169*t-409)*(409*t-169));
H2,phi:=ReducedMinimalWeierstrassModel(H);
J:=Jacobian(H2);
pJ:=Points(J : Bound:=20000);
inf:=[Q : Q in pJ | Order(Q) eq 0];
print "NPTS_J:", #pJ, "INFINITE_ORDER:", #inf, "t:", Cputime(T0);
if #inf gt 0 then
  Q:=inf[1];
  pts2:=Chabauty(Q);
  pts:=[Inverse(phi)(R) : R in pts2];
  print "CHABAUTY_NPTS:", #pts; print "PTS:", pts; print "t:", Cputime(T0);
end if;
