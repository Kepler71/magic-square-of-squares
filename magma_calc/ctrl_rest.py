from chab_all_lib import magma
cases=[("79/110",79,31,189),("79/110",79,110,189),("11/142",131,142,153),("11/142",142,11,153),("11/142",153,11,142),("48/163",115,48,163)]
for sl,c,a,b in cases:
    code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+{c}*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2));
J:=Jacobian(C); rb:=RankBound(J); print "RANKBOUND:", rb;
pts:=Points(C : Bound:=5000); print "SEARCH:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];
if rb eq 0 then print "CH0:", #Chabauty0(J); end if;
if rb eq 1 then
  cand:=[R : R in pts | R[3] ne 0 and R[2] ne 0];
  inf:=[R : R in pts | R[3] eq 0];
  Q:=[J!(R-inf[1]) : R in cand] cat [J!(R-S) : R in pts, S in pts | R ne S];
  Q:=[q : q in Q | Order(q) eq 0];
  if #Q gt 0 then S:=Chabauty(Q[1]); print "CHAB:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; else print "CHAB: нет точки бесконечного порядка среди разностей"; end if;
end if;"""
    print("===",sl,c,a,b); print(magma(code),flush=True)
