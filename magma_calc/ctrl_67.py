from chab_all_lib import magma
code="""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+235*z)*(1-67^2*z^2)*(1-101^2*z^2));
J:=Jacobian(C); rb:=RankBound(J); print "RB:", rb;
pts:=Points(C : Bound:=2000); print "SEARCH:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];
if rb eq 1 then Q:=[q : q in [J!(R-S) : R in pts, S in pts | R ne S] | Order(q) eq 0]; S:=Chabauty(Q[1]); print "CHAB:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; end if;"""
print(magma(code))
