# Проверка CTP: (1) симметрия, нулевая диагональ, кручение в ядре, независимость от локальных точек;
# (2) кривые E0/Q с полным 2-кручением, расширенные до k = Q(sqrt d): rank E(k) = r(E0) + r(E0^d) (PARI, точно),
#     избыток = dim Sel^2(E/k) - 2 - rank; должно быть rank CTP <= избыток (и = избыток, если Ш[2^oo] элементарна).
import sys, time
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else int(5))

def qrank(E):
    r = pari(E.ainvs()).ellinit().ellrank(); return ZZ(r[0]), ZZ(r[1])

def sel_elements(Cv, Sel):
    out = []
    for w in Sel.basis():
        a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n)); b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
        out.append((Cv.k(a), Cv.k(b)))
    return out

def ctp_matrix(Cv, els, reps=2):
    ct = CTP(Cv); ct.randomize_conic = True; ct.randomize_conic = True
    n = len(els); M = matrix(GF(2), n, n)
    for i in range(n):
        vals, nP = ct.pair(els[i], els, reps=reps)
        for j in range(n): M[i, j] = vals[j]
    return M

def run_case(d, rts, reps=2):
    k = QuadraticField(d, 'w')
    E0 = EllipticCurve(QQ, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)])
    r0 = qrank(E0); rd = qrank(E0.quadratic_twist(d))
    if r0[0] != r0[1] or rd[0] != rd[1]: return None
    rk = r0[0] + rd[0]
    Cv = Curve3(k, rts)
    Sel = Cv.selmer_full()
    exc = Sel.dimension() - 2 - rk
    if exc < 2: return ('skip', d, rts, Sel.dimension(), rk)
    els = sel_elements(Cv, Sel)
    t0 = time.time()
    M = ctp_matrix(Cv, els, reps=reps)
    sym = (M == M.transpose()); diag0 = all(M[i, i] == 0 for i in range(M.nrows()))
    # кручение: образы T1, T2 должны лежать в ядре
    e1, e2, e3 = Cv.e
    tors = [Cv.kappa(e1), Cv.kappa(e2)]
    ct = CTP(Cv); ct.randomize_conic = True
    # глобальные точки: генераторы E0(Q) и твиста E0^(d)(Q) (PARI ellrank) -> E(k)
    gl = []
    for P_ in pari(E0.ainvs()).ellinit().ellrank()[3]:
        x_ = QQ(P_[0]); gl.append(Cv.kappa(k(x_)))
    a1, a2, a3, a4, a6 = E0.ainvs()
    Et = EllipticCurve(QQ, [0, a2*d, 0, a4*d^2, a6*d^3])       # d y^2 = f(x):  X = d x
    for P_ in pari(Et.ainvs()).ellinit().ellrank()[3]:
        x_ = QQ(P_[0]) / d; gl.append(Cv.kappa(k(x_)))
    pts_ok = all(all(v == 0 for v in ct.pair(el, tors + gl, reps=1)[0]) for el in els) if gl else None
    return ('case', d, rts, Sel.dimension(), rk, exc, M.rank(), sym, diag0, pts_ok, len(gl), round(time.time() - t0))

if __name__ == '__main__' or True:
    cases = []
    for d in [2, 3, 5, 6, 7, 10, 13]:
        for a in range(1, 40):
            for b in range(1, 40):
                if gcd(a, b) != 1 or a == b: continue
                cases.append((d, [0, a, -b]))
    random.shuffle(cases)
    found = 0
    for d, rts in cases:
        try:
            res = run_case(d, rts)
        except Exception as ex:
            print("ERR", d, rts, type(ex).__name__, str(ex)[:120]); continue
        if res is None: continue
        if res[0] == "skip": print("skip", res[1:]); continue
        print(res); found += 1
        if found >= int(sys.argv[2] if len(sys.argv) > 2 else 8): break
