import sys, time
load('/home/kep/magicKube/descent/ek_descent.sage')
random.seed(int(1))
def qrank(E):
    r = pari(E.ainvs()).ellinit().ellrank(); return (ZZ(r[0]), ZZ(r[1]))
tests = [(5, [0, 1, -3]), (13, [0, 2, -7]), (15, [0, 5, -11]), (6, [-1, 0, 4])]
for d, rts in tests:
    k = QuadraticField(d, 'w')
    E0 = EllipticCurve(QQ, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)])
    r0 = qrank(E0); rd = qrank(E0.quadratic_twist(d))
    t0 = time.time()
    C = Curve3(k, rts)
    Sel = C.selmer_full()
    iso = [C.isogeny_descent(i) for i in range(3)]
    print(f"k=Q(sqrt{d}) E0 roots {rts}: rank E0(Q) in {r0}, twist in {rd} => rank E(k) in [{r0[0]+rd[0]},{r0[1]+rd[1]}];"
          f" full Sel2 dim {Sel.dimension()} -> bound {Sel.dimension()-2}; isogeny (s1,s2) {iso} -> bounds {[s+t-2 for s,t in iso]}  {time.time()-t0:.0f}s")
