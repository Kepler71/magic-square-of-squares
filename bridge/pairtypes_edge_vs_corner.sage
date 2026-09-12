print("6. sections: edge-fixed (proved) vs corner-fixed (untouched)")
secs=[(17,7,13),(7,1,5),(23,7,17),(31,17,25),(41,1,29),(47,23,37),(49,31,41),(73,17,53),(71,49,61),(89,23,65),(79,47,65)]
def sqf(x):
    x=ZZ(x); return sign(x)*prod(p^(e%2) for p,e in x.abs().factor())
print("  section        A_e    C_e   k_edge      | A_c    C_c    k_corner     AC_c sq? AC_e sq?")
for (b0,h0,n0) in secs:
    A=ZZ((h0^2+n0^2)/2); C=ZZ((b0^2+n0^2)/2); dk=sqf((C-A)*(C+A))
    Ac=ZZ(h0^2); Cc=ZZ(b0^2); dkc=sqf((Cc-Ac)*(Cc+Ac))
    print("  (%2d,%2d,%2d)  %6d %6d  Q(sqrt %4d) | %5d %6d  Q(sqrt %5d)   %5s   %5s"%(b0,h0,n0,A,C,dk,Ac,Cc,dkc,(Ac*Cc).is_square(),(A*C).is_square()))
print("  identity  k_edge = Q(sqrt(b0^2-h0^2)),  k_corner = Q(sqrt(2*(b0^2-h0^2))) for all 11:",
      all(sqf(((b0^2-h0^2)/2)*(2*n0^2))==sqf(b0^2-h0^2) and sqf((b0^2-h0^2)*(b0^2+h0^2))==sqf(2*(b0^2-h0^2)) for (b0,h0,n0) in secs))

print()
print("7. the elimination step X=(Y3-Y4)/(Y1-Y2): does the shift coefficient B enter?")
var('A C B p q Y1 Y2 Y3 Y4 Xx')
# Y1^2=A*p+B*q, Y2^2=A*p-B*q, Y3^2=C*p+B*q, Y4^2=C*p-B*q  (p=(t^2+1)^2 a square, q=t(t^2-1))
S=Y1+Y2; D=Y1-Y2
# X = (Y3-Y4)/(Y1-Y2) = (Y1+Y2)/(Y3+Y4) since (Y1-Y2)(Y1+Y2)=2Bq=(Y3-Y4)(Y3+Y4)
# u=S^2, v=D^2 :  u+v=4Ap ;  u/X^2+X^2*v=4Cp   =>
u = 4*p*Xx^2*(A*Xx^2-C)/(Xx^4-1); v = 4*p*(C*Xx^2-A)/(Xx^4-1)
print("   solve u+v=4Ap, u/X^2+X^2 v=4Cp :")
print("     u+v-4Ap        simplifies to", (u+v-4*A*p).simplify_full())
print("     u/X^2+X^2v-4Cp simplifies to", (u/Xx^2+Xx^2*v-4*C*p).simplify_full())
print("   => S^2 = 4p X^2 (A X^2 - C)/(X^4-1),  D^2 = 4p (C X^2 - A)/(X^4-1);  B does NOT appear.")
print("   so Q: (A X^2-C)(X^4-1)=sq  and  Q': (C X^2-A)(X^4-1)=sq  depend only on (A,C).")
print("   cross-check with VERIFY_CODEX_CHAIN_s15 (A=109,C=229): D^2 = 4(t^2+1)^2(229X^2-109)/(X^4-1) ->",
      (v.subs(A=109,C=229)).simplify_full())
