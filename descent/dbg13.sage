load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/run_sections.sage') if False else None
random.seed(int(3))
k = QuadraticField(30, 'r'); r = k.gen()
rts = [13906, 2312*r, -13906]
Cv = Curve3(k, rts)
P = [pl for pl in Cv.places if isinstance(pl, LocalSq) and pl.P == k.ideal(13, r + 2)][0]
a, b = Cv.isogeny_data(1); a2, b2 = -2*a, a^2 - 4*b
print("a,b =", a, b, " b2 = a^2-4b =", b2, " sqrt:", b2.is_square())
X = Cv.Rx.gen()
for name, g, bb in [("E", X*(X^2 + a*X + b), b), ("E'", X*(X^2 + a2*X + b2), b2)]:
    rt = [t for t, _ in g.roots()]
    print(name, "roots", rt, "valuations at P", [t.valuation(P.P) if t != 0 else 'inf' for t in rt], "coords", [P.coords(bb if t == 0 else t) for t in rt])
    rows = []; found = 0
    for it in range(20000):
        j = random.randint(-8, 10)
        x = random.choice(rt + [0]) + rand_elt(k, 13^6) * P.pi^j
        gx = g(x)
        if gx == 0 or not P.is_sq(gx): continue
        found += 1; rows.append(P.coords(x))
    print("  local points found", found, "image rank", matrix(GF(2), rows).rank() if rows else 0, "distinct", set(tuple(r_) for r_ in rows))
print("disc valuations:", [ (e1 - e2).valuation(P.P) for e1 in map(k, rts) for e2 in map(k, rts) if e1 != e2])
