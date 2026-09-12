import ast
rows = [ast.literal_eval(l) for l in open('/home/kep/magicKube/descent/ctp_cases_all.txt')]
seen = set()
for r in rows:
    d, rts = r[1], r[2]
    key = (d, tuple(rts))
    if key in seen: continue
    seen.add(key)
    E0 = EllipticCurve(QQ, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)])
    Et = E0.quadratic_twist(d)
    sh = []
    for E in (E0, Et):
        E = E.minimal_model()
        r_ = E.rank()
        s_ = E.sha().an() if r_ <= 1 else None
        sh.append((r_, s_))
    two = [ (s_ // 2^0).valuation(2) if s_ else None for (_, s_) in sh ]
    print(f"d={d} roots={rts}: k-excess {r[5]}, CTP rank {r[6]} | E0: rank {sh[0][0]}, Sha_an {sh[0][1]} | twist: rank {sh[1][0]}, Sha_an {sh[1][1]}"
          f" | v2(|Sha_Q|) total = {sum(t for t in two if t is not None)}")
