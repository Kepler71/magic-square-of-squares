# Claude, 26.09.2026. Дополнительная проверка своей 4-клеточной модели: для наклона найти 4-клеточные T (T ≠ −T,
# по одному из пары T, −T) ранга 1 и прогнать check_factor. Это второе, независимое по модели закрытие.
import sys, os, json, time, itertools
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rigdem
from rigdem import *

rigdem.LOGF = open(os.path.join(here, 'ctrl_4cell.log'), 'a')
slopes = sys.argv[1].split(',')
NUSE = int(sys.argv[2]) if len(sys.argv) > 2 else 2
TLIM = float(sys.argv[3]) if len(sys.argv) > 3 else 500
res_all = {}
t0 = time.time()
for sl in slopes:
    r, s = map(int, sl.split('/'))
    S = sorted(cells(r, s))
    cand = []
    for T in itertools.combinations(S, 4):
        T = list(T)
        negT = sorted(-x for x in T)
        if negT == T or negT < T:
            continue
        cand.append(T)
    log('%s: 4-клеточных T (по модулю ±, T ≠ −T): %d' % (sl, len(cand)))
    used = []
    ranks = []
    for T in cand:
        if time.time() - t0 > TLIM or len(used) >= NUSE:
            break
        md = model(T)
        if md['E1'].discriminant() == 0:
            continue
        er = pari(md['E1'].minimal_model()).ellrank()
        rk = [int(er[0]), int(er[1])]
        ranks.append([T, rk])
        if rk != [1, 1]:
            continue
        res = check_factor(r, s, T, ctrl_brute=400)
        D = md['E1'].discriminant()
        ps = [p for p in prime_range(1009, 3000) if D % p != 0 and all(
            QQ(c).denominator() % p != 0 for tt in md['E1'].torsion_points() if not tt.is_zero() for c in tt.xy())][:2]
        res['algebra_ctrl_Fp'] = {str(p): algebra_ctrl_Fp(T, p) for p in ps}
        short = {k: res.get(k) for k in ('ellrank', 'mwrank_rank', 'tors_order', 'tors_gcd_Fp', 'index_bound', 'own_psat_all', 'pari_sat_same',
                                          'model_used', 'B_rig', 'B_silverman_sage', 'c1', 'C', 'hP0_lower', 'hP0_sage', 'norm_ok', 'M0',
                                          'nondegenerate', 'caseB_zero_poly', 'closed_algebraic', 'aux_ok', 'closed')}
        log('  T=%s: %s' % (T, json.dumps(short, ensure_ascii=False)))
        log('    перебор: %s; за M0: %s; F_p: %s' % (res['ctrl_brute'], res['ctrl_beyond'], json.dumps(res['algebra_ctrl_Fp'])))
        used.append(res)
    res_all[sl] = dict(ranks=ranks, results=used)
json.dump(res_all, open(os.path.join(here, 'ctrl_4cell.json'), 'w'), ensure_ascii=False, indent=1, default=str)
log('ctrl_4cell: готово за %.0f с' % (time.time() - t0))
