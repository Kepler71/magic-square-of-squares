# Claude, 26.09.2026 (сводка STATUS). Своя быстрая перепроверка единственного «трудного» наклона вне 45: 96/407.
# Для каждого 3-клеточного T ⊂ Λ = {±s, ±r, ±(s−r), ±(s+r)} (с точностью до T ~ −T):
#   y² = ∏(1+λz)  ⇔  Y² = ∏(X + L/λ), X = Lz, Y = Ly, L = ∏λ.
# Если PARI ellrank даёт верхнюю границу 0, то E(Q) = кручение; все z = X/L проверяются на «восемь клеток — квадраты».
# Читал до расчёта: запись fable_bl/extra_checks.json['96/407'] (там 5 множителей ранга 0) и build_master.py::nondeg.
from sage.all import *
import itertools, time
r, s = 96, 407
LAM = sorted({r, s - r, s, s + r} | {-r, -(s - r), -s, -(s + r)})
def is_sq(q):
    q = QQ(q); return q >= 0 and q.numerator().is_square() and q.denominator().is_square()
def nondeg(z):
    return z != 0 and all(is_sq(1 + l * z) for l in LAM)
seen = set(); t0 = time.time(); closing = []; failed = []
for T in itertools.combinations(LAM, 3):
    key = min(tuple(sorted(T)), tuple(sorted(-x for x in T)))
    if key in seen: continue
    seen.add(key)
    L = prod(T)
    R = PolynomialRing(QQ, 'X'); X = R.gen()
    f = prod(X + QQ(L) / l for l in T)
    if f.discriminant() == 0: continue
    E = EllipticCurve([0, f[2], 0, f[1], f[0]])
    E = E.global_minimal_model()
    # ellrank через отдельный процесс gp (в cypari2 на этом наклоне segfault, как и в census_300)
    import subprocess
    a = [int(c) for c in E.a_invariants()]
    pr = subprocess.run(['gp', '-q', '-f', '-D', 'parisizemax=2000000000'], input=f'print(ellrank(ellinit({a}))[1..2])\n',
                         capture_output=True, text=True, timeout=300)
    out = pr.stdout.strip()
    try:
        lo, hi = [int(v) for v in out.strip('[]').split(',')]
    except ValueError:
        print(f'  T={T} a={a}: gp не дал ответа (rc={pr.returncode}): {pr.stderr.strip()[:200]!r}', flush=True); failed.append(T); continue
    if hi != 0: continue
    tors = E.torsion_points()
    Ew = EllipticCurve([0, f[2], 0, f[1], f[0]])
    tors = Ew.torsion_points()
    zs = [('inf' if P.is_zero() else P[0] / L) for P in tors]
    bad = [z for z in zs if z != 'inf' and nondeg(z)]
    closing.append((T, int(lo), int(hi), len(tors), bad))
    print(f'[{time.time()-t0:.1f}s] T={T} ellrank=[{lo},{hi}] #tors={len(tors)} z={zs} невырожденных={bad}', flush=True)
print('классов T:', len(seen), 'ранга 0 (верх. граница PARI):', len(closing),
      'из них без невырожденных z:', sum(1 for c in closing if not c[4]))
print('gp не посчитал:', failed)
print('96/407 закрыт этим методом:', any(not c[4] for c in closing))
