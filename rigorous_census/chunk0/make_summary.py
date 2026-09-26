# Claude, 26.09.2026. Сводка по итоговым rig_*.json + сравнение с joint_sieve/dem_<r>_<s>.json (Fable).
# Пишет summary_chunk0.json и table_chunk0.md. Ничего не пересчитывает, кроме сравнения чисел.
import json, glob

base = '/home/kep/magicKube/rigorous_census/chunk0'
SLOPES = '143/206 96/211 61/217 19/221 211/230 204/247 231/250 65/261 43/278 265/298 108/301 251/308 237/317 140/319 197/325'.split()
key = lambda T: min(tuple(sorted(T)), tuple(sorted(-t for t in T)))

ver = json.load(open(f'{base}/verify_chunk0.json'))
sat = {(v['slope'], tuple(v['T'])): v for v in json.load(open(f'{base}/sat10k.json'))}
summ = {}; rows = []; allc = []
for sl in SLOPES:
    r, s = sl.split('/')
    d = json.load(open(f'{base}/rig_{r}_{s}.json'))
    fab = json.load(open(f'/home/kep/magicKube/joint_sieve/dem_{r}_{s}.json'))
    res = d['results']
    closers = [o for o in res if o.get('closed_rigorous')]
    allc += [(sl, o) for o in closers]
    vset = {tuple(v['T']): v for v in ver[sl]['details']}
    for o in closers:
        assert vset[tuple(o['T'])]['ok'] and vset[tuple(o['T'])]['M0'] == o['M0_rig']
        assert sat[(sl, tuple(o['T']))]['sat_p_le_1e4']
        assert o['ctrl_ineq_fail'] == 0 and o['ctrl_T_square_fail'] == 0 and o['norm_ok']
        assert o['caseA_nondeg'] == [] and o['caseB_nondeg'] == []
    mw1 = [o for o in closers if o.get('mwrank_rank_bound') == 1]
    best = min(mw1, key=lambda o: (o['M0_rig'], -o['hP0_lo']))
    mine = {key(o['T']): o for o in res}
    fcmp = []
    for fo in fab['results']:
        o = mine[key(fo['T'])]
        fcmp.append(dict(T=fo['T'], closed=bool(o.get('closed_rigorous')), mwrank=o.get('mwrank_rank_bound'),
                         M0_rig=o.get('M0_rig'), M0_fable=fo['M0'],
                         hP0_fable=fo['hP0'], hP0_in_bracket=bool(o['hP0_lo'] - 1e-9 <= fo['hP0'] <= o['hP0_hi'] + 1e-9),
                         B_fable=fo['B'], B_S_mine=o['B_silverman'], B_rig=o['B_rig'],
                         c1_fable=fo['c1'], c1_mine=o['c1'], fable_ok=fo['ok'], fable_nondeg=fo['nondeg_solutions']))
    summ[sl] = dict(
        n_classes=len(res), n_rank11=sum(o.get('ellrank') == [1, 1] for o in res),
        n_rank11_nopoint=sum('без явной' in str(o.get('status')) for o in res),
        n_closed=len(closers), n_closed_mwrank1=len(mw1),
        n_errors=sum(o.get('status') == 'ОШИБКА' for o in res),
        best=dict(T=best['T'], hP0=[best['hP0_lo'], best['hP0_hi']], B_rig=best['B_rig'], B_S=best['B_silverman'],
                  c1=best['c1'], C=best['C_rig'], M0=best['M0_rig'], M_run=best['M_run'], cand=best['caseA_cand'],
                  tors=best['tors_struct'], caseB_roots_z=best['caseB_roots_z']),
        fable=fcmp, closed_rigorous=bool(closers), nondeg_solutions=[])
    b = best
    fstr = ', '.join(('✓' if f['closed'] else '✗') + f"{f['M0_rig']}/{f['M0_fable']}" for f in fcmp)
    rows.append(f"| {sl} | {len(closers)} ({len(mw1)}) | {b['T']} | {b['hP0_lo']:.3f} | {b['B_rig']:.2f} / {b['B_silverman']:.2f} "
                f"| {b['c1']:.2f} | {b['M0_rig']} | {b['caseA_cand']} | 0 | {fstr} |")
ratios = [o['B_rig'] / o['B_silverman'] for _, o in allc]
g = dict(
    n_slopes_closed=sum(v['closed_rigorous'] for v in summ.values()),
    n_classes=sum(v['n_classes'] for v in summ.values()),
    n_rank11=sum(v['n_rank11'] for v in summ.values()),
    n_rank11_nopoint=sum(v['n_rank11_nopoint'] for v in summ.values()),
    n_closers=len(allc), n_closers_mwrank1=sum(v['n_closed_mwrank1'] for v in summ.values()),
    min_closers_mwrank1_per_slope=min(v['n_closed_mwrank1'] for v in summ.values()),
    n_errors=sum(v['n_errors'] for v in summ.values()),
    cand_total=sum(o['caseA_cand'] for _, o in allc),
    M0_range=[min(o['M0_rig'] for _, o in allc), max(o['M0_rig'] for _, o in allc)],
    B_rig_range=[min(o['B_rig'] for _, o in allc), max(o['B_rig'] for _, o in allc)],
    B_up_range=[min(o['B_up'] for _, o in allc), max(o['B_up'] for _, o in allc)],
    B_S_range=[min(o['B_silverman'] for _, o in allc), max(o['B_silverman'] for _, o in allc)],
    ratio_Brig_BS=[min(ratios), max(ratios)], all_Brig_gt_BS=all(o['Brig_gt_BS'] for _, o in allc),
    all_Brig_from_K=all(o['B_low'] >= o['B_up'] for _, o in allc),
    bracket_width_max=max(o['hP0_hi'] - o['hP0_lo'] for _, o in allc),
    all_M0rig_le_M0fablestyle=all(o['M0_rig'] <= o['M0_fable_style'] for _, o in allc),
    fable_45_closed=sum(f['closed'] for v in summ.values() for f in v['fable']),
    fable_M0rig_le_M0fable=sum(f['M0_rig'] is not None and f['M0_rig'] <= f['M0_fable'] for v in summ.values() for f in v['fable']),
    fable_hP0_in_bracket=sum(f['hP0_in_bracket'] for v in summ.values() for f in v['fable']),
    fable_B_equals_my_BS=sum(abs(f['B_fable'] - f['B_S_mine']) < 1e-6 for v in summ.values() for f in v['fable']),
    fable_c1_less_than_mine=[(sl, f['T'], round(f['c1_fable'], 3), round(f['c1_mine'], 3)) for sl, v in summ.items()
                             for f in v['fable'] if f['c1_fable'] < f['c1_mine'] - 1e-9],
    caseB_nonempty=sum(bool(o['caseB_roots_z']) for _, o in allc),
    tors_structs=sorted(set(str(o['tors_struct']) for _, o in allc)),
)
json.dump(dict(global_=g, slopes=summ), open(f'{base}/summary_chunk0.json', 'w'), ensure_ascii=False, indent=1)
hdr = ('| наклон | закрывающих (mwrank = 1) | лучший T | ĥ(P₀) ≥ | B строг. / B_S | c₁ | M₀ | кандидатов | z ≠ 0 | Fable: M₀ строг./его |\n'
       '|---|---|---|---:|---|---:|---:|---:|---:|---|')
open(f'{base}/table_chunk0.md', 'w').write(hdr + '\n' + '\n'.join(rows) + '\n')
print(hdr); print('\n'.join(rows)); print(json.dumps(g, ensure_ascii=False, indent=1))
