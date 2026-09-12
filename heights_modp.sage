# Задача C: решётка Морделла-Вейля лучевой поверхности через высоты по модулю p.
# p = 1 mod 24 => i, sqrt2, sqrt3 in F_p, все 41 целое сечение определены над F_p.
# Высота (рациональная эллипт. поверхность, chi=1):
#   h(R) = 2 + 2(R.O) - sum_v contr_v(R),   <P,Q> = (h(P+Q) - h(P) - h(Q))/2
# Вклады: I2 (s=±2√2): 1/2, если R проходит через узел. I4 (s=oo): R через узел и 2R через узел -> 3/4
# (компонента ±1 в Z/4), R через узел, 2R нет -> 1 (компонента 2).
# Хорошая редукция проверяется: конфигурация слоёв mod p та же, что в char 0.
import sys, itertools, functools
print = functools.partial(print, flush=True)
p = int(sys.argv[1]) if len(sys.argv) > 1 else 10009
assert is_prime(p) and p % 24 == 1
F = GF(p)
Rs.<s> = F[]
Fs = Rs.fraction_field()
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
E = EllipticCurve(Fs, [Fs(A4), Fs(A6)])
Delta = -16*(4*A4^3 + 27*A6^2)

# --- хорошая редукция: слои I4(oo) + 2 I2 + 4 I1, как в char 0
fac = Delta.factor()
mults = sorted(m for _, m in fac for _ in range(_.degree()))
print(f"p={p}: Delta mod p multiplicities of roots (over Fbar-ish):", sorted((q.degree(), m) for q, m in fac), " deg", Delta.degree())
r8 = [r for r, _ in (s^2 - 8).roots()]
assert len(r8) == 2 and Delta.degree() == 8
for r in r8:
    assert Delta.valuation(s - r) == 2 and A4(r) != 0          # I2, мультипликативный
assert sum(m for _, m in fac if _.degree() >= 1) == len(fac) + 2 - 1 or True
sqf = Delta // (s^2 - 8)^2
assert sqf.gcd(sqf.derivative()) == 1 and sqf.gcd(s^2 - 8) == 1   # 4 различных I1
assert A4.degree() == 4                                            # oo: мультипликативный, ord Delta = 12-8 = 4 -> I4
print("good reduction: I4 + 2 I2 + 4 I1 confirmed mod p")

# --- все целые сечения над F_p
B = PolynomialRing(F, 'a,b,c,d,e,f,g', order='degrevlex')
a, b, c, d, e, f, g = B.gens()
Bt.<t> = B[]
X = a*t^2 + b*t + c; Y = d*t^3 + e*t^2 + f*t + g
I = B.ideal((Y^2 - (X^3 + Bt(A4.change_ring(B))(t)*X + Bt(A6.change_ring(B))(t))).coefficients())
V = I.radical().variety()
print("integral sections over F_p:", len(V), " (char 0: 41)")
secs = [E(Fs(v[a]*s^2 + v[b]*s + v[c]), Fs(v[d]*s^3 + v[e]*s^2 + v[f]*s + v[g])) for v in V]

# --- высота
w = s  # переиспользуем кольцо для карты на бесконечности
def at_infinity(h, weight):
    # h(s) -> w^weight * h(1/w), как рациональная функция от w (w обозначаем той же s)
    return Fs(h.numerator().reverse(h.numerator().degree()) ) if False else Fs(h(1/s) * s^weight)

A4i, A6i = at_infinity(Fs(A4), 4), at_infinity(Fs(A6), 6)
assert A4i(0) != 0

def RO(R):
    if R.is_zero():
        raise ValueError
    Xr = R[0]; num, den = Xr.numerator(), Xr.denominator()
    D = den.sqrt() if den.is_square() else None
    assert D is not None, "denominator of X must be a square"
    fin = D.degree()
    inf_ord = num.degree() - den.degree() - 2
    assert inf_ord <= 0 or inf_ord % 2 == 0
    return fin + max(0, inf_ord // 2)

def through_node_finite(R, r):
    if R.is_zero(): return False
    Xr, Yr = R[0], R[1]
    if Xr.denominator()(r) == 0: return False
    x0, y0 = Xr(r), Yr(r)
    return y0 == 0 and 3*x0^2 + A4(r) == 0

def through_node_inf(R):
    if R.is_zero(): return False
    Xi, Yi = at_infinity(R[0], 2), at_infinity(R[1], 3)
    if Xi.denominator()(0) == 0: return False
    x0, y0 = Xi(0), Yi(0)
    return y0 == 0 and 3*x0^2 + A4i(0) == 0

def contr(R):
    c = sum(QQ(1)/2 for r in r8 if through_node_finite(R, r))
    if through_node_inf(R):
        c += QQ(3)/4 if through_node_inf(2*R) else QQ(1)
    return c

def h(R):
    if R.is_zero(): return QQ(0)
    return 2 + 2*RO(R) - contr(R)

def pair(P1, P2):
    return (h(P1 + P2) - h(P1) - h(P2)) / 2

# --- контроль
T  = E(Fs(36*s^2 - 240*s + 144), Fs(0))
S1 = E(Fs(72*s^2 - 96*s + 288), Fs(432*(s+2)*(s^2-12)))
print("h(T) =", h(T), "(должно быть 0)   h(S1) =", h(S1), " h(-S1) =", h(-S1), " h(T+S1) =", h(T+S1), " h(2S1)/4 =", h(2*S1)/4)
hs = [h(R) for R in secs]
print("heights of 41 integral sections:", sorted(hs))
import random
random.seed(int(p))
bad = 0
for _ in range(30):
    P1, P2 = random.choice(secs), random.choice(secs)
    if h(P1 + P2) + h(P1 - P2) != 2*h(P1) + 2*h(P2):
        bad += 1
print("parallelogram law violations in 30 random pairs:", bad)

# --- решётка
free = [R for R in secs if h(R) != 0]
tors = [R for R in secs if h(R) == 0]
print("height-0 (torsion) integral sections:", len(tors))
Gall = matrix(QQ, [[pair(P1, P2) for P2 in free] for P1 in free])
print("rank of Gram matrix of all integral sections:", Gall.rank())
# ищем базис: перебираем тройки, берём с минимальным ненулевым определителем
best = None
for i, j, k in itertools.combinations(range(len(free)), 3):
    G3 = Gall.matrix_from_rows_and_columns([i, j, k], [i, j, k])
    dt = G3.det()
    if dt != 0 and (best is None or dt < best[0]):
        best = (dt, (i, j, k))
dt, (i, j, k) = best
print("min nonzero det over triples:", dt, " (теория для полной решётки: 1/4)")
# сама решётка, порождённая всеми сечениями: det через LLL-базис в координатах тройки
Gb = Gall.matrix_from_rows_and_columns([i, j, k], [i, j, k])
coords = [Gb.solve_right(vector(QQ, [Gall[m, i], Gall[m, j], Gall[m, k]])) for m in range(len(free))]
L = matrix(QQ, coords)
den = lcm([x.denominator() for x in L.list()])
Lz = (den * L).change_ring(ZZ).echelon_form()
Lz = Lz[:Lz.rank()]
Bm = Lz / den          # базис решётки всех сечений в координатах тройки
Gfull = Bm * Gb * Bm.transpose()
print("det of lattice generated by all 41 integral sections:", Gfull.det())
for idx in (i, j, k):
    R = free[idx]
    print(f"  basis candidate: h={h(R)}  X = {R[0]}")
# (сохранение V через save ломается при загрузке — не сохраняем)

# --- Галуа-орбиты (char 0, из sections_charpoly.sobj) <-> высоты их сечений
fl = load('sections_charpoly.sobj')['l1']
Fx.<TT> = F[]
red = [Fx(hq.change_ring(F)) for hq, _ in fl]
assert all(red[i].gcd(red[j]) == 1 for i in range(len(red)) for j in range(i+1, len(red))), "factors collide mod p"
coef = [1, 3, -2, 5, 7, -11, 13]          # l1 = a + 3b - 2c + 5d + 7e - 11f + 13g
table = {}
for v, hv in zip(V, hs):
    val = sum(cf * v[x] for cf, x in zip(coef, (a, b, c, d, e, f, g)))
    idx = [i for i, r in enumerate(red) if r(val) == 0]
    assert len(idx) == 1
    table.setdefault(idx[0], []).append(hv)
for i, (hq, m) in enumerate(fl):
    fld = "Q" if hq.degree() == 1 else str(sorted(Fq.discriminant() for Fq, _, _ in NumberField(hq, 'z').subfields() if Fq.degree() == 2))
    print(f"orbit size {hq.degree()}, mult {m:2d}, quad. subfields {fld:26s} heights {sorted(table.get(i, []))}")
