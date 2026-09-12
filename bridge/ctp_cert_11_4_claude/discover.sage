# -*- coding: utf-8 -*-
# Шаг 1: квартика g1 целевого класса delta=(1,274,274); квартики для базиса Sel^2;
#        спаривание Касселса-Тейта <g1, g2> по Теореме 3.1 Фишера.
# Ранг E(Q) НЕ используется нигде.
import sys, json, time
load('/home/kep/magicKube/bridge/ctp_cert_11_4_claude/ctp_quartic_snapshot.sage')

D = '/home/kep/magicKube/bridge/ctp_cert_11_4_claude/'
(m, n, s, b, e, c, rootsM, I, J, phis, order, delta_trip) = load(D + 'setup.sobj')
R.<X> = QQ[]
fM = prod(X - r for r in rootsM)
M = EllipticCurve([0, 0, 0, fM.coefficients(sparse=False)[1], fM.coefficients(sparse=False)[0]])

F = FisherCTP(QQ, I, J, verbose=True)
def mkdelta(trip):
    return F.E.from_comps([QQ(trip[i]) for i in order])

# ------------------------------------------------ группа Селмера (список из bridge/selmer_11_4.log)
SEL = [(-28770,2,-14385),(-28770,137,-210),(-2055,210,-1918),(-2055,14385,-7),
       (-959,1,-959),(-959,274,-14),(-274,105,-28770),(-274,28770,-105),
       (-105,2,-210),(-105,137,-14385),(-30,210,-7),(-30,14385,-1918),
       (-14,1,-14),(-14,274,-959),(-1,105,-105),(-1,28770,-28770),
       (1,1,1),(1,274,274),(14,105,30),(14,28770,2055),
       (30,2,15),(30,137,4110),(105,210,2),(105,14385,137),
       (274,1,274),(274,274,1),(959,105,2055),(959,28770,30),
       (2055,2,4110),(2055,137,15),(28770,210,137),(28770,14385,2)]
PR = [-1, 2, 3, 5, 7, 137]
def vec(trip):
    v = []
    for d in trip:
        d = QQ(d); u = d.numerator()*d.denominator()
        row = [1 if u < 0 else 0]
        u = abs(u)
        for p in PR[1:]:
            k2 = 0
            while u % p == 0:
                u //= p; k2 += 1
            row.append(k2 % 2)
        assert u == 1, (trip, u)
        v += row
    return vector(GF(2), v)

VS = VectorSpace(GF(2), 18)
W = VS.subspace([vec(t) for t in SEL])
print("dim Sel^2 (из списка) =", W.dimension(), " |Sel^2| =", 2^W.dimension())
assert W.dimension() == 5
TARGET = (1, 274, 274)
assert vec(TARGET) in W

# базис Sel^2, содержащий целевой класс
basis = [TARGET]
cur = VS.subspace([vec(TARGET)])
for t in SEL:
    if cur.dimension() == 5:
        break
    if vec(t) not in cur:
        basis.append(t); cur = VS.subspace([vec(x) for x in basis])
print("базис Sel^2:", basis)

# --------------------------------------------- произведение классов
def mul(t1, t2):
    out = []
    for a, b2 in zip(t1, t2):
        p = QQ(a)*QQ(b2)
        # приведение по модулю квадратов
        num = p.numerator()*p.denominator()
        sgn = -1 if num < 0 else 1
        num = abs(num); red = sgn
        for p2, k2 in factor(num):
            if k2 % 2:
                red *= p2
        out.append(QQ(red))
    return tuple(out)

# --------------------------------------------- построение квартик с кэшем
CACHE = D + 'quartics.sobj'
try:
    Qcache = load(CACHE)
except Exception:
    Qcache = {}

def quartic(trip):
    key = str(tuple(trip))
    if key in Qcache:
        return [QQ(x) for x in Qcache[key]]
    t0 = time.time()
    g = F.quartic_from_delta(mkdelta(trip), tag=key)
    print("  квартика для %s построена за %.1f c: %s" % (key, time.time()-t0, g))
    sys.stdout.flush()
    Qcache[key] = [str(x) for x in g]
    save(Qcache, CACHE)
    return g

print("\n=== квартика целевого класса ===")
g1 = quartic(TARGET)
print("g1 =", g1)
print("I(g1), J(g1) =", quartic_I(g1), quartic_J(g1))
assert quartic_I(g1) == I and quartic_J(g1) == J
z1 = F.z_inv(g1)
print("z(g1) =", z1, " компоненты:", [F.E.comp(z1, i) for i in range(3)])
zc = [F.E.comp(z1, i) for i in range(3)]
# класс z(g1) в порядке (e1,e2,e3)
cls = [None]*3
for i in range(3):
    cls[order[i]] = zc[i]
def sqfree(a):
    a = QQ(a); num = a.numerator()*a.denominator()
    sgn = -1 if num < 0 else 1; num = abs(num); red = sgn
    for p2, k2 in factor(num):
        if k2 % 2:
            red *= p2
    return red
print("класс z(g1) в (Q*/Q*^2)^3, порядок (e1,e2,e3):", [sqfree(x) for x in cls])
assert tuple(sqfree(x) for x in cls) == TARGET, "z(g1) НЕ в целевом классе!"
print(">>> z(g1) лежит ровно в требуемом классе (1,274,274)")
sys.stdout.flush()

# --------------------------------------------- спаривание с базисом Sel^2
print("\n=== спаривание <g1, h> по базису Sel^2 ===")
results = {}
for t2 in basis:
    if tuple(t2) == TARGET:
        pass
    t3 = mul(TARGET, t2)
    print("\n-- g2 класс %s, g3 класс %s" % (str(t2), str(t3)))
    sys.stdout.flush()
    try:
        g2 = quartic(t2)
        g3 = quartic(t3)
    except Exception as ex:
        print("   ОШИБКА построения:", repr(ex)); sys.stdout.flush(); continue
    if g2[0] == 0 and g3[0] == 0:
        print("   g2(1,0)=g3(1,0)=0 — тривиальный класс, спаривание 0"); continue
    t0 = time.time()
    try:
        val, npl = F.pair(g1, g2, g3, reps=2)
    except Exception as ex:
        print("   ОШИБКА спаривания:", repr(ex)); sys.stdout.flush(); continue
    print("   <g1,g2> = %s  (мест %d, %.1f c)" % (val, npl, time.time()-t0))
    results[str(tuple(t2))] = int(val)
    sys.stdout.flush()
    if val:
        print("   >>> НЕНУЛЕВОЕ СПАРИВАНИЕ НАЙДЕНО")
        save((g1, g2, g3, tuple(t2), tuple(t3)), D + 'found.sobj')

print("\nИТОГ:", results)
save(results, D + 'results.sobj')
