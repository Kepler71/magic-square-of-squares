import json
R.<X> = QQ[]
res = {}
for lbl, ai in [("1088e2",[0,0,0,-364,-2640]), ("582d2",[1,0,0,-194,-1056]), ("1025b2",[1,-1,0,-667,2616])]:
    E = EllipticCurve(ai)
    Es = E.short_weierstrass_model()
    a4 = Es.a4(); a6 = Es.a6()
    rts = sorted([r[0] for r in (X^3 + a4*X + a6).roots(QQ)])
    if len(rts) < 3:
        print(lbl, "не полное 2-кручение в short model"); continue
    I = -48*a4; J = -1728*a6
    cov = pari(Es).ell2cover()
    qs = []
    for c in cov:
        co = [QQ(t) for t in pari(c[0]).Vec()]
        while len(co) < 5: co = [QQ(0)] + co
        qs.append([str(t) for t in co])
    res[lbl] = {"ai": [str(t) for t in ai], "short": [str(a4), str(a6)],
                "roots": [str(t) for t in rts], "I": str(I), "J": str(J),
                "rank": int(E.rank()), "sha": float(E.sha().an_numerical(prec=50)),
                "quartics": qs}
    print(lbl, "a4,a6 =", a4, a6, " roots", rts, " I,J", I, J, " nq", len(qs))
    for q in qs: print("   ", q)
json.dump(res, open("ctrl_curves.json", "w"), indent=1)
print("ok")
