d = 3; rts = [0, 19, -10]
E0 = EllipticCurve(QQ, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)]).minimal_model()
Et = E0.quadratic_twist(d).minimal_model()
for name, E in [("E0", E0), ("E0^(3)", Et)]:
    r = pari(E.ainvs()).ellinit().ellrank()
    ar = E.analytic_rank()
    s = E.sha().an() if ar <= 1 else None
    print(name, E.ainvs(), "cond", E.conductor(), "PARI ellrank", r[:3], "analytic rank", ar, "Sha_an", s,
          "2-Selmer rank", E.selmer_rank(), "torsion", E.torsion_order())
