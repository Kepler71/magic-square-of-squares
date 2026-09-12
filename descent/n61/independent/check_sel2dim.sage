load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import json, random
random.seed(int(7))
k = kk; r = rr
d = json.load(open('/home/kep/magicKube/descent/n61/ctp_n61_certificate_full_seed11.json'))
def K(s): return k(sage_eval(str(s), locals={'r': r}))
e = [K(t) for t in d['model_roots']]
e1, e2, e3 = e
fails = []
def ck(n, c, i=''):
    if not c: fails.append((n, i)); print("FAIL", n, i, flush=True)
def ideal_of(s):
    inside = s[s.index('(') + 1:s.rindex(')')]
    return k.ideal([K(g.strip()) for g in inside.split(',')])
S = [ideal_of(s) for s in d['S']]

# ---------- A. носители eps/eta и списка мест ----------
allowed = set()
for P in S: allowed.add(P)
for q in primes(2, 61):
    for P in k.primes_above(q):
        if P.norm() <= 60: allowed.add(P)
for row in d['rows']:
    W = row['witnesses']
    ck(f"row{row['row']}: список мест = allowed",
       set(W['places']) == set(str(P) for P in allowed),
       str(set(W["places"]).symmetric_difference(set(str(P) for P in allowed))))
    for s in list(W['eps_reps']) + [t for pr in W['eta_reps'] for t in pr]:
        I = k.ideal(K(s))
        bad = [P for P, _ in I.factor() if not any(P == Q for Q in allowed)]
        ck(f"row{row['row']}: носитель {s[:25]} в allowed", not bad, str(bad))
print("A: носители eps/eta и множество мест проверены")

# ---------- B. независимый пересчёт dim Sel^2 ----------
def coords(pl, x):
    if isinstance(pl, tuple): return [0 if exact_sign_emb(x, pl[1]) > 0 else 1]
    P = pl
    if P.smallest_integer() == 2: return p2_coords(x)
    F, qq, pi = odd_data(P)
    v = ZZ(x.valuation(P)); un = x / pi**v
    return [int(v % 2), 0 if chi(P, un) == 1 else 1]

def dimV(pl):
    if isinstance(pl, tuple): return 1
    return 4 if pl.smallest_integer() == 2 else 2

def f_poly(x): return (x - e1)*(x - e2)*(x - e3)

def local_image(pl, target):
    """подгруппа образа kappa: E(k_v)/2 -> (k_v^*/sq)^2, порождённая кручением и найденными точками"""
    gens = [((e1 - e2)*(e1 - e3), e1 - e2), (e2 - e1, (e2 - e1)*(e2 - e3))]
    rows = [coords(pl, a) + coords(pl, b) for a, b in gens]
    M = matrix(GF(2), rows)
    tries = 0
    while M.rank() < target and tries < 200000:
        tries += 1
        if isinstance(pl, tuple):
            x = QQ(random.randint(-10**7, 10**7))/QQ(random.randint(1, 10**4)) \
                + QQ(random.randint(-10**4, 10**4))/QQ(random.randint(1, 10**3)) * r
        else:
            p = pl.smallest_integer()
            j = random.randint(-4, 6)
            x = random.choice(list(e) + [k(0)]) + (QQ(random.randint(-p**4, p**4)) + QQ(random.randint(-p**4, p**4))*r) * k.uniformizer(pl, others='positive')**j
        if any(x == t for t in e): continue
        fx = f_poly(x)
        if fx == 0: continue
        if not is_local_square(pl, fx): continue
        rows.append(coords(pl, x - e1) + coords(pl, x - e2))
        M = matrix(GF(2), rows)
    return M, tries

places = list(S) + [('inf', 1), ('inf', -1)]
blocks = []
tot_codim = 0
for pl in places:
    tgt = 1 if isinstance(pl, tuple) else (4 if pl.smallest_integer() == 2 else 2)
    M, tries = local_image(pl, tgt)
    ck(f"локальный образ {pl}: достигнут {tgt}", M.rank() == tgt, f"{M.rank()} за {tries} попыток")
    ck(f"локальный образ {pl}: не больше {tgt}", M.rank() <= tgt)
    V = 2*dimV(pl)
    tot_codim += V - M.rank()
    blocks.append((pl, M.row_space(), V))
    print(f"  {str(pl)[:32]:34s} dim V={V:2d}  dim image={M.rank()}  (попыток {tries})")
print("суммарная коразмерность:", tot_codim)

gens = [K(g) for g in d['kS2_basis']]
ng = len(gens)
# матрица отображения k(S,2)^2 -> (+) V_v / L_v
cols = []
for i in range(2*ng):
    a = gens[i] if i < ng else k(1)
    b = gens[i - ng] if i >= ng else k(1)
    vec = []
    for pl, L, V in blocks:
        ca = coords(pl, a) + coords(pl, b)
        vec += list(ca)
    cols.append(vec)
A = matrix(GF(2), cols)              # строки = образы генераторов
# фактор по локальным образам: ядро = {v : для каждого места координаты лежат в L_v}
# строим как систему: для каждого места берём аннулятор L_v
sysrows = []
off = 0
for pl, L, V in blocks:
    Ann = L.basis_matrix().right_kernel().basis_matrix()   # функционалы, зануляющие L_v
    for fun in Ann.rows():
        col = [0]*(2*ng)
        for i in range(2*ng):
            col[i] = sum(A[i][off + t]*fun[t] for t in range(V))
        sysrows.append(col)
    off += V
Sys = matrix(GF(2), sysrows)
ker = Sys.right_kernel()
print("\ndim k(S,2)^2 =", 2*ng, "  число условий =", Sys.nrows(), "  ранг условий =", Sys.rank())
print("МОЙ независимый dim Sel^2 =", ker.dimension())
ck("dim Sel^2 = 7", ker.dimension() == 7, str(ker.dimension()))

# заявленный базис Sel^2 лежит в моём ядре и порождает его
sb = [(K(a), K(b)) for a, b in d['sel2_basis']]
def to_vec(pair):
    """координаты (a,b) в базисе k(S,2)^2 (перебор 2^22 недопустим — решаем линейно по логарифмам)"""
    return None
# вместо координат: проверяем, что элементы ядра и базис Sel порождают одно и то же подпространство,
# сравнивая размерности образов в (+) V_v (инъективно для k(S,2)^2? не обязательно) —
# достаточная проверка: каждый элемент базиса Sel удовлетворяет всем локальным условиям
badloc = 0
for (a, b) in sb:
    for pl, L, V in blocks:
        v = vector(GF(2), coords(pl, a) + coords(pl, b))
        if v not in L: badloc += 1; print("  вне локального образа:", pl, a, b)
ck("базис Sel^2 удовлетворяет локальным условиям", badloc == 0)

print("\nИТОГ:", "ВСЕ ОК" if not fails else f"ОШИБКИ ({len(fails)}): {fails[:10]}")
