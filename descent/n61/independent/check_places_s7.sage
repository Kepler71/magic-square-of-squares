load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import json, sys, time
k = kk; r = rr
CERT = '/home/kep/magicKube/descent/n61/ctp_n61_certificate_seed7.json'
d = json.load(open(CERT))
def K(s): return k(sage_eval(str(s), locals={'r': r}))
e = [K(t) for t in d['model_roots']]
R4 = PolynomialRing(k, 'u0,u1,u2,u3'); u = R4.gens()
fails = []
def ck(n, c, i=''):
    if not c:
        fails.append((n, i))
        print("FAIL", n, i, flush=True)

# ---- интервальные знаки L_j в вещественном месте (собственная реализация) ----
def real_signs(vals, Lc, sr):
    """vals = [(x-e_i)/eps_i] > 0; Lc[j] = список коэффициентов; точка (1, sqrt v1, sqrt v2, sqrt v3)."""
    for prec in [256, 1024, 4096, 16384, 65536, 262144]:
        RI = RealIntervalField(prec); sD = RI(165).sqrt() * sr
        def ev(a):
            a0, a1 = list(k(a)); return RI(a0) + RI(a1) * sD
        us = [RI(1)] + [ev(v_).sqrt() for v_ in vals]
        out = {}
        ok = True
        for j in (0, 1, 2):
            S = sum(ev(Lc[j][i]) * us[i] for i in range(4))
            if S.contains_zero(): ok = False; break
            out[j] = 1 if S > 0 else -1
        if ok: return out, prec
    raise RuntimeError("знак не отделён")

t0 = time.time()
Mrows = []
for row in d['rows']:
    i = row['row']; W = row['witnesses']
    eps = [K(t) for t in W['eps_reps']]
    eta = [(K(a), K(b)) for a, b in W['eta_reps']]
    Lc = {j: [K(W['L_forms'][str(j)][str(t)]) for t in range(4)] for j in range(4)}
    recs = row['places']
    nrec = len(recs)
    # разбиение на повторы: (27 конечных + 2 вещественных)
    fin_places = W['places']
    blk = len(fin_places) + 2
    ck(f"row{i}: число записей", nrec % blk == 0, f"{nrec} / {blk}")
    nrep = nrec // blk
    totals = []
    for rep in range(nrep):
        tot = [0]*len(eta)
        chunk = recs[rep*blk:(rep+1)*blk]
        # конечные места
        for ridx, rc in enumerate(chunk):
            if str(rc['place']).startswith('real'):
                sr = 1 if str(rc['place']).endswith('=1') else -1
                x = K(rc['x'])
                vals = [K(v) for v in rc['vals_under_sqrt']]
                # точность: vals должны быть ровно (x - e_i)/eps_i и положительны
                for t in range(3):
                    ck(f"row{i}/rep{rep}/real{sr}: val{t}", vals[t] == (x - e[t])/eps[t])
                    ck(f"row{i}/rep{rep}/real{sr}: val{t}>0", exact_sign_emb(vals[t], sr) > 0)
                sg, prec = real_signs(vals, Lc, sr)
                rec_sg = {j: int(rc['signs_L'][str(j)]) for j in (0, 1, 2)}
                ck(f"row{i}/rep{rep}/real{sr}: знаки L", sg == rec_sg, f"{sg} vs {rec_sg}")
                f1s = sg[1]*sg[0]; f2s = sg[2]*sg[0]
                for t, (h1, h2) in enumerate(eta):
                    s1 = -1 if (f1s < 0 and exact_sign_emb(h2, sr) < 0) else 1
                    s2 = -1 if (f2s < 0 and exact_sign_emb(h1, sr) < 0) else 1
                    rs1, rs2 = int(rc['symbols'][t][0]), int(rc['symbols'][t][1])
                    ck(f"row{i}/rep{rep}/real{sr}/eta{t}: s1", s1 == rs1, f"{s1} vs {rs1}")
                    ck(f"row{i}/rep{rep}/real{sr}/eta{t}: s2", s2 == rs2, f"{s2} vs {rs2}")
                    tot[t] += 0 if s1*s2 == 1 else 1
                continue
            # ---- конечное место ----
            P = k.ideal([K(g) for g in rc['gens_of_place']])
            ck(f"row{i}/rep{rep}/rec{ridx}: место совпадает со списком", str(P) == fin_places[ridx],
               f"{P} vs {fin_places[ridx]}")
            p = P.smallest_integer()
            e2 = ZZ(k(2).valuation(P))
            N = ZZ(rc['N'])
            x = K(rc['x'])
            us = [K(t) for t in rc['u_approx']]
            ck(f"row{i}/rep{rep}/{P}: u0=1", us[0] == 1)
            A = [rc['A_i'][t] for t in range(4)]
            ck(f"row{i}/rep{rep}/{P}: A0=inf", str(A[0]) == '+Infinity')
            Amine = [Infinity]
            for t in range(3):
                a_t = (x - e[t]) / eps[t]
                ck(f"row{i}/rep{rep}/{P}: a{t}!=0", a_t != 0)
                va = ZZ(a_t.valuation(P))
                ck(f"row{i}/rep{rep}/{P}: v(a{t}) чётна", va % 2 == 0, str(va))
                # заявленная относительная точность корня
                dv = (us[t+1]**2 - a_t).valuation(P)
                ck(f"row{i}/rep{rep}/{P}: u{t+1}^2 ~ a{t}", dv >= N + va, f"v(u^2-a)={dv} < N+v(a)={N+va}")
                # Гензель: a_t — локальный квадрат
                ck(f"row{i}/rep{rep}/{P}: a{t} лок.квадрат (Гензель)", dv - va > 2*e2, f"{dv-va} <= {2*e2}")
                Ai = N - e2 + va/2
                Amine.append(Ai)
                ck(f"row{i}/rep{rep}/{P}: A{t+1}", ZZ(A[t+1]) == Ai, f"{A[t+1]} vs {Ai}")
            # значения L
            Lv = {}
            for j in (0, 1, 2):
                val = sum(Lc[j][t]*us[t] for t in range(4))
                ck(f"row{i}/rep{rep}/{P}: L{j} значение", val == K(rc['L_values'][str(j)]))
                ck(f"row{i}/rep{rep}/{P}: L{j}!=0", val != 0)
                Lv[j] = val
                vL = ZZ(val.valuation(P))
                ck(f"row{i}/rep{rep}/{P}: v(L{j})", vL == ZZ(rc['v_L'][str(j)]), f"{vL} vs {rc['v_L'][str(j)]}")
                AL = min(ZZ(Lc[j][t].valuation(P)) + Amine[t] for t in range(1, 4) if Lc[j][t] != 0)
                marg = AL - vL - 2*e2
                ck(f"row{i}/rep{rep}/{P}: margin{j}", marg == ZZ(rc['margins'][str(j)]),
                   f"{marg} vs {rc['margins'][str(j)]}")
                ck(f"row{i}/rep{rep}/{P}: margin{j}>0", marg > 0, str(marg))
            f1 = Lv[1]/Lv[0]; f2 = Lv[2]/Lv[0]
            ck(f"row{i}/rep{rep}/{P}: f1", f1 == K(rc['f1']))
            ck(f"row{i}/rep{rep}/{P}: f2", f2 == K(rc['f2']))
            for t, (h1, h2) in enumerate(eta):
                s1 = hsym(P, f1, h2); s2 = hsym(P, f2, h1)
                rs1, rs2 = int(rc["symbols"][t][0]), int(rc["symbols"][t][1])
                ck(f"row{i}/rep{rep}/{P}/eta{t}: s1", s1 == rs1, f"mine={s1} cert={rs1}  f1={f1} h2={h2}")
                ck(f"row{i}/rep{rep}/{P}/eta{t}: s2", s2 == rs2, f"mine={s2} cert={rs2}  f2={f2} h1={h1}")
                tot[t] += 0 if s1*s2 == 1 else 1
        totals.append([t % 2 for t in tot])
    ck(f"row{i}: повторы согласованы", all(t == totals[0] for t in totals), str(totals))
    rec_vals = [int(v) for v in row['values']]
    ck(f"row{i}: значения строки", totals[0] == rec_vals, f"{totals[0]} vs {rec_vals}")
    Mrows.append(totals[0])
    print(f"row {i}: {totals[0]}   ({time.time()-t0:.0f}s)", flush=True)

n = len(d['rows'])
M = matrix(GF(2), [row[:n] for row in Mrows])
print("\nМоя матрица CTP:"); print(M)
Mc = matrix(GF(2), d['matrix'])
ck("матрица совпадает с сертификатом", M == Mc)
ck("симметрия", M == M.transpose())
ck("нулевая диагональ", all(M[i, i] == 0 for i in range(n)))
print("ранг:", M.rank())
print("столбцы T1,T2,G (должны быть нулевыми):", [[row[n+j] for row in Mrows] for j in range(3)])
ck("T1,T2,G в ядре", all(row[n+j] == 0 for row in Mrows for j in range(3)))
print("граница ранга: rank E(k) <=", n - 2 - M.rank())
print("\nИТОГ:", "ВСЕ ОК" if not fails else f"ОШИБКИ ({len(fails)})")
if fails:
    for f_ in fails[:40]: print("  ", f_)
