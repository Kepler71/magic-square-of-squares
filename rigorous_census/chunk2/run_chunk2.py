# Claude, 26.09.2026. Прогон строгой проверки по 15 наклонам куска 2 (последовательно, один процесс).
# Множители T — те же три, что у Fable в joint_sieve/dem_<r>_<s>.json (берётся только список T;
# ранг, образующая, кручение, границы и перебор — свои). Если ни один не закрыл — перебор своих T.
import sys, os, json, time, itertools
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rigdem
from rigdem import *

SLOPES = ['165/439', '422/441', '425/441', '333/442', '111/445', '329/449', '438/449', '126/451',
          '172/451', '295/452', '305/461', '236/475', '165/493', '52/499', '129/499']
BRUTE = int(sys.argv[1]) if len(sys.argv) > 1 else 300
only = sys.argv[2].split(',') if len(sys.argv) > 2 else None

rigdem.LOGF = open(os.path.join(here, 'run_chunk2.log'), 'a')
log('=== старт: %d наклонов, контроль перебором |a|,b ≤ %d' % (len(SLOPES), BRUTE))
summary = {}
for sl in SLOPES:
    if only and sl not in only:
        continue
    r, s = map(int, sl.split('/'))
    fab = json.load(open('/home/kep/magicKube/joint_sieve/dem_%d_%d.json' % (r, s)))
    Tlist = [x['T'] for x in fab['results']]
    log('наклон %s: множители %s' % (sl, Tlist))
    results = []
    for T in Tlist:
        try:
            res = check_factor(r, s, T, ctrl_brute=BRUTE)
        except Exception:
            import traceback
            res = dict(T=T, error=traceback.format_exc()[-2000:])
        # контроль алгебры случая B над F_p
        if 'error' not in res and 'skip' not in res:
            md = model(T); D = md['E1'].discriminant()
            ps = [p for p in prime_range(1009, 3000) if D % p != 0 and all(
                QQ(c).denominator() % p != 0 for tt in md['E1'].torsion_points() if not tt.is_zero() for c in tt.xy())][:2]
            res['algebra_ctrl_Fp'] = {str(p): algebra_ctrl_Fp(T, p) for p in ps}
        short = {k: res.get(k) for k in ('ellrank', 'mwrank_rank', 'sat_index_eclib', 'pari_sat_same', 'tors_order', 'tors_gcd_Fp',
                                          'model_used', 'B_rig', 'B_silverman_sage', 'c1', 'C', 'hP0_lower', 'hP0_sage', 'norm_ok', 'M0',
                                          'nA', 'n_values', 'index_bound', 'own_psat_all', 'nondegenerate', 'caseB_zero_poly', 'closed_algebraic', 'aux_ok', 'closed', 'skip', 'error', 'sec')}
        log('  T=%s: %s' % (T, json.dumps(short, ensure_ascii=False, default=str)))
        if 'ctrl_brute' in res:
            log('    контроль перебором: симметричных %s, пропущено %s' % (res['ctrl_brute']['symmetric_found'], res['ctrl_brute']['missing']))
        if 'algebra_ctrl_Fp' in res:
            log('    контроль за M0: %s' % json.dumps(res.get('ctrl_beyond'), ensure_ascii=False))
            log('    контроль алгебры B над F_p: %s' % json.dumps(res['algebra_ctrl_Fp']))
        results.append(res)
    closed = [x['T'] for x in results if x.get('closed')]
    alert = [x['nondegenerate'] for x in results if x.get('nondegenerate')]
    json.dump(dict(slope=sl, results=results, closed_by=closed, closed_rigorous=bool(closed), alert=alert),
              open(os.path.join(here, 'rig_%d_%d.json' % (r, s)), 'w'), ensure_ascii=False, indent=1, default=str)
    summary[sl] = dict(closed_by=closed, alert=alert)
    log('наклон %s: закрыт строго множителями %s; ALERT %s' % (sl, closed, alert))
json.dump(summary, open(os.path.join(here, 'summary%s.json' % ('' if not only else '_part')), 'w'), ensure_ascii=False, indent=1)
log('=== конец')
