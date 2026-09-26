# Claude, 26.09.2026. Сводка куска 2 по двум проходам: rig_*.json (Sage-проход, rigdem.py) и indep/indep_*.json
# (чисто питоновский проход, indep/indep_verify.py). Наклон «закрыт строго», если есть множитель T, закрытый в ОБОИХ
# проходах с одинаковым M0, с рангом ≤ 1 по PARI ellrank, и (для «двух реализаций ранга») eclib rank_bound = 1.
import json, os
here = os.path.dirname(os.path.abspath(__file__))
SL = ['165/439', '422/441', '425/441', '333/442', '111/445', '329/449', '438/449', '126/451',
      '172/451', '295/452', '305/461', '236/475', '165/493', '52/499', '129/499']
out = {}
for sl in SL:
    r, s = map(int, sl.split('/'))
    A = {tuple(x['T']): x for x in json.load(open(os.path.join(here, 'rig_%d_%d.json' % (r, s))))['results']}
    Bj = {tuple(x['T']): x for x in json.load(open(os.path.join(here, 'indep', 'indep_%d_%d.json' % (r, s))))['results']}
    fac = []
    for T in A:
        a = A[T]; b = Bj.get(T, {})
        row = dict(T=list(T), closed_pass1=bool(a.get('closed')), closed_pass2=bool(b.get('closed')),
                   M0_pass1=a.get('M0'), M0_pass2=b.get('M0'), model=a.get('model_used'),
                   B_rigorous=round(a['B_rig'], 4), B_silverman_sage=round(a['B_silverman_sage'], 4),
                   B_rig_gt_silverman=a['B_rig'] > a['B_silverman_sage'],
                   c1=round(a['c1'], 4), C=round(a['C'], 4), hP0_lower=round(min(a['hP0_lower'], b.get('hP0_lo', 1e9)), 5),
                   hP0_sage=round(a['hP0_sage'], 5), ellrank=a['ellrank'], mwrank_rank_bound_pass2=b.get('mwrank_bound'),
                   tors=a['tors_order'], index_bound=a['index_bound'], sat_primes_own=sorted(b.get('saturation_own', {}).keys(), key=int),
                   nondegenerate=a.get('nondegenerate', []) + b.get('nondegenerate', []),
                   brute_bound_pass1=a.get('ctrl_brute', {}).get('bound'), brute_missing=a.get('ctrl_brute', {}).get('missing', []) + b.get('ctrl_brute', {}).get('missing', []))
        row['both'] = row['closed_pass1'] and row['closed_pass2'] and row['M0_pass1'] == row['M0_pass2'] and not row['nondegenerate']
        row['two_rank_impl'] = row['ellrank'] == [1, 1] and row['mwrank_rank_bound_pass2'] == 1
        fac.append(row)
    closed_by = [f['T'] for f in fac if f['both']]
    out[sl] = dict(status='закрыт строго (ПО)' if closed_by else 'не удалось', closed_by=closed_by,
                   closed_by_with_two_rank_impl=[f['T'] for f in fac if f['both'] and f['two_rank_impl']], factors=fac)
    print(sl, out[sl]['status'], 'множителей %d/3, из них с двумя реализациями ранга %d' % (len(closed_by), len(out[sl]['closed_by_with_two_rank_impl'])),
          'M0:', [f['M0_pass1'] for f in fac], 'B:', [f['B_rigorous'] for f in fac], 'Bsil:', [f['B_silverman_sage'] for f in fac])
json.dump(out, open(os.path.join(here, 'final_chunk2.json'), 'w'), ensure_ascii=False, indent=1)
print('закрыто строго: %d из %d' % (sum(1 for v in out.values() if v['closed_by']), len(out)))
