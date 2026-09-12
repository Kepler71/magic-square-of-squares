load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import json
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
print("S =", [ (P.smallest_integer(), P.norm()) for P in S])

# --- 1. S содержит все плохие простые: 2 и делители (e_i - e_j) ---
need = set()
for P in k.primes_above(2): need.add(P)
for (a, b) in [(0,1),(0,2),(1,2)]:
    for P, _ in k.ideal(e[a] - e[b]).factor(): need.add(P)
missing = [P for P in need if not any(P == Q for Q in S)]
ck("S покрывает 2 и (e_i-e_j)", not missing, str(missing))
print("плохие простые:", sorted([(P.smallest_integer(), P.norm()) for P in need]))

# --- 2. базис k(S,2) ---
gens = [K(g) for g in d['kS2_basis']]
print("dim k(S,2) заявлен:", len(gens))
# независимость по модулю квадратов
bad = 0
for mask in range(1, 2**len(gens)):
    pass  # 2^11 = 2048 — считаем ниже через линейную алгебру валентностей+знаков
# ожидаемая размерность: |S| + rk O^* + 1 = 9 + 1 + 1 = 11
ck("dim k(S,2) = |S|+rk(O*)+1", len(gens) == len(S) + 1 + 1, str(len(gens)))
# независимость: никакое нетривиальное произведение не квадрат (перебор 2^11)
sq = 0
for mask in range(1, 2**len(gens)):
    p = k(1)
    for i in range(len(gens)):
        if mask >> i & 1: p *= gens[i]
    if p.is_square(): sq += 1; print("  зависимость, mask", mask); break
ck("базис k(S,2) независим", sq == 0)

# --- 3. образы кручения и точки G ---
t1 = ((e1 - e2)*(e1 - e3), e1 - e2)
t2 = (e2 - e1, (e2 - e1)*(e2 - e3))
rec_t = [(K(a), K(b)) for a, b in d['torsion_images']]
ck("kappa(T1)", (rec_t[0][0]/t1[0]).is_square() and (rec_t[0][1]/t1[1]).is_square(), str(rec_t[0]))
ck("kappa(T2)", (rec_t[1][0]/t2[0]).is_square() and (rec_t[1][1]/t2[1]).is_square(), str(rec_t[1]))
q = 175015061936
Eo = EllipticCurve(k, [0, -9747472064*r, 0, -q^2, 9747472064*r*q^2])
G = Eo(-173017629136, -784866450593280*r - 4316765478263040)
xG = G[0]/244^2
fG = (xG - e1)*(xG - e2)*(xG - e3)
ck("G на кривой (f(xG) квадрат)", fG.is_square(), str(fG))
ck("y(G)/244^3 согласовано", (G[1]/244^3)^2 == fG)
recG = (K(d['G_image'][0]), K(d['G_image'][1]))
ck("kappa(G)", ((xG - e1)/recG[0]).is_square() and ((xG - e2)/recG[1]).is_square(), str(recG))
# независимость T1,T2,G в k^*/(k^*)^2 x k^*/(k^*)^2
trip = [rec_t[0], rec_t[1], recG]
dep = 0
for mask in range(1, 8):
    a = k(1); b = k(1)
    for i in range(3):
        if mask >> i & 1: a *= trip[i][0]; b *= trip[i][1]
    if a.is_square() and b.is_square(): dep += 1; print("  зависимость T/G, mask", mask)
ck("T1,T2,G независимы", dep == 0)

# --- 4. базис Sel^2 независим ---
sb = [(K(a), K(b)) for a, b in d['sel2_basis']]
dep = 0
for mask in range(1, 2**len(sb)):
    a = k(1); b = k(1)
    for i in range(len(sb)):
        if mask >> i & 1: a *= sb[i][0]; b *= sb[i][1]
    if a.is_square() and b.is_square(): dep += 1; print("  зависимость Sel, mask", mask); break
ck("базис Sel^2 независим (dim 7)", dep == 0)
# T1,T2,G лежат в <Sel basis>?
def in_span(pair):
    for mask in range(2**len(sb)):
        a = k(1); b = k(1)
        for i in range(len(sb)):
            if mask >> i & 1: a *= sb[i][0]; b *= sb[i][1]
        if (a/pair[0]).is_square() and (b/pair[1]).is_square(): return True
    return False
for nm, pr in [('T1', rec_t[0]), ('T2', rec_t[1]), ('G', recG)]:
    ck(f"{nm} в Sel^2", in_span(pr))

# --- 5. локальные образы 2-спуска: собственные координаты и ранги ---
def coords(place, x):
    if isinstance(place, tuple):
        return [0 if exact_sign_emb(x, place[1]) > 0 else 1]
    P = place
    if P.smallest_integer() == 2:
        return p2_coords(x)
    F, qq, pi = odd_data(P)
    v = ZZ(x.valuation(P))
    un = x / pi**v
    return [int(v % 2), 0 if chi(P, un) == 1 else 1]

print("\nместо                         target  achieved(cert)  achieved(мой)  теория")
for w in d['sel2_local_witness']:
    nm = w['place']
    if nm == 'real':
        pls = [('inf', 1), ('inf', -1)]
        theory = 1
    else:
        P = ideal_of(nm); pls = [P]
        theory = 2 if P.smallest_integer() != 2 else 2 + ZZ(P.ramification_index()*P.residue_class_degree())
    reps = [(K(a), K(b)) for a, b in w['image_reps']]
    for pl in pls:
        rows = [coords(pl, a) + coords(pl, b) for a, b in reps]
        rk = matrix(GF(2), rows).rank() if rows else 0
        mark = 'OK' if (rk == int(w['target_dim']) == theory == int(w['achieved_dim'])) else '<<<'
        print(f"{str(pl)[:28]:30s} {w['target_dim']:>5}  {w['achieved_dim']:>13}  {rk:>12}  {theory:>6}  {mark}")
        ck(f"локальный образ {pl}: достигнут target", rk == int(w['target_dim']))
        ck(f"локальный образ {pl}: target = теория", int(w['target_dim']) == theory)
    ck(f"cert achieved == target для {nm}", int(w['achieved_dim']) == int(w['target_dim']))

print("\nИТОГ:", "ВСЕ ОК" if not fails else f"ОШИБКИ ({len(fails)}): {fails[:10]}")
