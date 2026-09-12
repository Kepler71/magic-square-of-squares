load('/home/kep/magicKube/descent/ek_descent.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
for name, A, C, D in [('s15', 109, 229, 65), ('n17', 169, 409, 34), ('n61', 3061, 4381, 671), ('n79', 3217, 5233, 455), ('n89', 2377, 6073, 910)]:
    k, rts = E1_roots(A, C, D)
    e = [k(t) for t in rts]
    out = []
    for i in range(3):
        ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
        b = (ei - ej)*(ei - ek_)
        out.append(bool(b.is_square()))
    E = EllipticCurve(k, [0, -sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], -prod(e)])
    print(f"{name}: b_i квадрат? {out}; изогенный класс: {len(E.isogeny_class().curves) if k.degree()==1 else '—'}")
