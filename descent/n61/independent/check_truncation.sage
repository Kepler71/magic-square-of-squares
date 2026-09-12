"""Для P вне списка из 27 мест: eps_i, eta_j и коэффициенты L_j — P-единицы, кривая C_eps имеет
хорошую редукцию. Ищем F_q-точку на C_eps с L_j != 0 для j=0,1,2; она поднимается по Гензелю,
значит f1, f2 — P-единицы, и (f1,eta2)_P = (f2,eta1)_P = 1. Проверяем это явно для N(P) <= BOUND."""
load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import json
k = kk; r = rr
d = json.load(open('/home/kep/magicKube/descent/n61/ctp_n61_certificate_full_seed11.json'))
def K(s): return k(sage_eval(str(s), locals={'r': r}))
e = [K(t) for t in d['model_roots']]
def ideal_of(s):
    inside = s[s.index('(') + 1:s.rindex(')')]
    return k.ideal([K(g.strip()) for g in inside.split(',')])
allowed = [ideal_of(s) for s in d['S']]
for q in primes(2, 61):
    for P in k.primes_above(q):
        if P.norm() <= 60 and not any(P == Q for Q in allowed): allowed.append(P)
BOUND = 4000
cand = []
for q in primes(2, BOUND + 1):
    for P in k.primes_above(q):
        if P.norm() <= BOUND and P.smallest_integer() != 2 and not any(P == Q for Q in allowed):
            cand.append(P)
print("простых вне списка с N(P) <=", BOUND, ":", len(cand))
fails = []
for row in d['rows']:
    i = row['row']; W = row['witnesses']
    eps = [K(t) for t in W['eps_reps']]
    Lc = {j: [K(W['L_forms'][str(j)][str(t)]) for t in range(4)] for j in range(4)}
    nbad = 0
    for P in cand:
        F = P.residue_field(); q = F.cardinality()
        eb = [F(t) for t in e]; epb = [F(t) for t in eps]
        if len(set(eb)) < 3 or any(t == 0 for t in epb): nbad += 1; fails.append((i, P, 'плохая редукция')); continue
        Lb = {j: [F(c) for c in Lc[j]] for j in range(3)}
        if any(all(c == 0 for c in Lb[j]) for j in range(3)):
            nbad += 1; fails.append((i, P, 'L_j == 0 mod P')); continue
        found = False
        for xb in F:
            vals = [(xb - eb[t]) / epb[t] for t in range(3)]
            if any(v == 0 for v in vals): continue
            if any(v**((q - 1) // 2) != F(1) for v in vals): continue
            rts = [v.sqrt() for v in vals]
            for s1 in (1, -1):
                for s2 in (1, -1):
                    for s3 in (1, -1):
                        pt = [F(1), s1*rts[0], s2*rts[1], s3*rts[2]]
                        if all(sum(Lb[j][t]*pt[t] for t in range(4)) != 0 for j in range(3)):
                            found = True; break
                    if found: break
                if found: break
            if found: break
        if not found: nbad += 1; fails.append((i, P, 'нет точки с L_j != 0'))
    print(f"row {i}: {len(cand) - nbad}/{len(cand)} простых подтверждено (символ = 1)", flush=True)
print("\nИТОГ:", "ВСЕ ОК" if not fails else f"ПРОБЛЕМЫ ({len(fails)}): {fails[:10]}")
