load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import json
k = kk; r = rr
CERT = '/home/kep/magicKube/descent/n61/ctp_n61_certificate_full_seed11.json'
d = json.load(open(CERT))
def K(s): return k(sage_eval(str(s), locals={'r': r}))

e = [K(t) for t in d['model_roots']]
e1, e2, e3 = e
print("roots:", e)
fails = []
def ck(n, c, i=''):
    if not c: fails.append((n, i)); print("FAIL", n, i)

R4 = PolynomialRing(k, 'u0,u1,u2,u3'); u = R4.gens()
S = [k.ideal(2)]  # заполним ниже реальным S
# множество "разрешённых" простых: S из сертификата + все простые нормы <= 60
allowed = set()
S_cert = []
for s in d['S']:
    # восстанавливаем идеал из строки вида "Fractional ideal (a, b)"
    inside = s[s.index('(') + 1:s.rindex(')')]
    gens = [K(g.strip()) for g in inside.split(',')]
    S_cert.append(k.ideal(gens))
for P in S_cert: allowed.add(P)
for q in primes(2, 61):
    for P in k.primes_above(q):
        if P.norm() <= 60: allowed.add(P)
print("|allowed primes| =", len(allowed))

sel_basis = [(K(a), K(b)) for a, b in d['sel2_basis']]
tors = [(K(a), K(b)) for a, b in d['torsion_images']]
Gim = (K(d['G_image'][0]), K(d['G_image'][1]))
targets = sel_basis + tors + [Gim]

for row in d['rows']:
    i = row['row']; W = row['witnesses']
    eps_rec = [K(t) for t in W['eps_reps']]
    ep1, ep2, ep3 = eps_rec
    # (a) eps_reps совпадает по классу с eps строки и eps3 = eps1*eps2
    e_row = [K(row['eps'][0]), K(row['eps'][1])]
    ck(f"row{i}: eps1 класс", (e_row[0] / ep1).is_square(), str(e_row[0] / ep1))
    ck(f"row{i}: eps2 класс", (e_row[1] / ep2).is_square(), str(e_row[1] / ep2))
    ck(f"row{i}: eps3=eps1eps2", ep3 == ep1 * ep2)
    # (b) eta_reps совпадают по классу с sel2_basis+torsion+G
    eta_rec = [(K(a), K(b)) for a, b in W['eta_reps']]
    ck(f"row{i}: len eta", len(eta_rec) == len(targets))
    for t, ((h1, h2), (g1, g2)) in enumerate(zip(eta_rec, targets)):
        ck(f"row{i}: eta{t} кл1", (h1 / g1).is_square(), f"{h1} / {g1}")
        ck(f"row{i}: eta{t} кл2", (h2 / g2).is_square(), f"{h2} / {g2}")
    # (c) коники из пучка
    G1 = ep1*u[1]**2 - ep2*u[2]**2 - (e2 - e1)*u[0]**2
    G2 = ep1*u[1]**2 - ep3*u[3]**2 - (e3 - e1)*u[0]**2
    cones = {3: G1, 2: G2, 1: G1 - G2, 0: (e3 - e1)*G1 - (e2 - e1)*G2}
    Lmine = {}
    for j, Q in cones.items():
        vars3 = [t for t in range(4) if t != j]
        ck(f"row{i}: конус{j} без u{j}", Q.monomial_coefficient(u[j]**2) == 0)
        ck(f"row{i}: конус{j} диагональна",
           all(Q.monomial_coefficient(u[a]*u[b]) == 0 for a in range(4) for b in range(a+1, 4)))
        cf = [Q.monomial_coefficient(u[t]**2) for t in vars3]
        cf_rec = [K(t) for t in W['conic_coeffs'][str(j)]]
        ck(f"row{i}: коэф коники {j}", cf == cf_rec, f"{cf} vs {cf_rec}")
        pt = [K(t) for t in W['conic_points'][str(j)]]
        val = sum(cf[t]*pt[t]**2 for t in range(3))
        ck(f"row{i}: точка на конике {j}", val == 0, str(val))
        ck(f"row{i}: точка {j} нетривиальна", any(t != 0 for t in pt))
        Lmine[j] = sum(2*cf[t]*pt[t]*u[vars3[t]] for t in range(3))
    # (d) сверка L с записанными с точностью до общего множителя
    for j in range(4):
        rec = {t: K(W['L_forms'][str(j)][str(t)]) for t in range(4)}
        mine = {t: Lmine[j].monomial_coefficient(u[t]) for t in range(4)}
        ck(f"row{i}: носитель L{j}", [t for t in range(4) if rec[t] != 0] == [t for t in range(4) if mine[t] != 0],
           f"{rec} vs {mine}")
        nz = [t for t in range(4) if mine[t] != 0]
        if nz:
            lam = rec[nz[0]] / mine[nz[0]]
            ck(f"row{i}: L{j} пропорц.", all(rec[t] == lam*mine[t] for t in nz), f"lam={lam}")
            # L_j должна занулять квадрику конуса вдоль образующей: проверка div — через то,
            # что L_j — касательная плоскость: подстановка p в L_j даёт 2*Q(p) = 0
            vars3 = [t for t in range(4) if t != j]
            pt = [K(t) for t in W['conic_points'][str(j)]]
            subs_pt = [k(0)]*4
            for tt in range(3): subs_pt[vars3[tt]] = pt[tt]
            ck(f"row{i}: L{j}(p)=0", Lmine[j](subs_pt) == 0, str(Lmine[j](subs_pt)))
        # (e) идеал коэффициентов нормированной формы поддержан в allowed
        cs = [rec[t] for t in range(4) if rec[t] != 0]
        I = k.ideal(cs)
        supp_ok = all(any(P == Q for Q in allowed) for P, _ in I.factor())
        ck(f"row{i}: носитель идеала L{j} в allowed", supp_ok, str(I.factor()))
    print(f"row {i}: коники/точки/формы проверены")

print("\nИТОГ:", "ВСЕ ОК" if not fails else f"ОШИБКИ ({len(fails)}): {fails[:10]}")
