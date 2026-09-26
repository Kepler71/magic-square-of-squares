# Claude, 26.09.2026: сводка ГЛАВНОЙ ТАБЛИЦЫ переписи наклонов 3 <= s <= 500 (без пересчёта кривых).
# Только чтение чужих файлов; пишет в rigorous_census/census_500_master.json и rigorous_census/master/.
# Все циклы ограничены: s <= 500, файлы конечны.
import json, re, os, sys, collections, glob, time
from math import gcd, isqrt
from fractions import Fraction as F
ROOT='/home/kep/magicKube'
OUT_JSON=f'{ROOT}/rigorous_census/census_500_master.json'
OUT_CHK=f'{ROOT}/rigorous_census/master/checks.json'
def P(*a): print(*a, flush=True)
def key(sl): r,s=map(int,sl.split('/')); return (s,r)
def rs(sl): r,s=map(int,sl.split('/')); return r,s
def Lam(r,s): v=[r,s-r,s,s+r]; return set(v)|{-x for x in v}
def issq(x):
    x=F(x); return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
def nondeg(z,r,s):
    # True, если z != 0 делает все восемь нецентральных клеток квадратами (т.е. НАСТОЯЩЕЕ решение)
    if z in ('inf',): return False
    t=F(z)
    if t==0: return False
    return all(issq(1+e*t) for e in Lam(r,s))
def jl(fn):
    out=[]
    for l in open(fn):
        l=l.strip()
        if l: out.append(json.loads(l))
    return out
t0=time.time()
checks=collections.OrderedDict(); problems=[]; notes=[]
def chk(name, ok, info=''):
    checks[name]={'ok':bool(ok),'info':info}
    P(('OK  ' if ok else 'FAIL'), name, info)
    if not ok: problems.append(f'{name}: {info}')

# ---------- 0. перечень наклонов
ALL=[f'{r}/{s}' for s in range(3,501) for r in range(1,s) if gcd(r,s)==1 and 2*r!=s]
RANGES=[('3-200',3,200),('201-300',201,300),('301-400',301,400),('401-500',401,500)]
def rng(sl):
    s=rs(sl)[1]
    for nm,a,b in RANGES:
        if a<=s<=b: return nm
cnt=collections.Counter(rng(x) for x in ALL)
chk('0.counts', len(ALL)==76114 and cnt['3-200']==12230 and cnt['201-300']==15166 and cnt['301-400']==21280 and cnt['401-500']==27438, f'{len(ALL)} {dict(cnt)}')
ALLSET=set(ALL)
ev=collections.defaultdict(list)   # slope -> list of evidence dicts
def add(sl, m, b, f, **kw):
    if sl not in ALLSET: problems.append(f'чужой наклон {sl} в {f}'); return
    d=dict(m=m,b=b,f=f); d.update(kw); ev[sl].append(d)

# ---------- 1. s <= 200
cov=json.load(open(f'{ROOT}/criterion_proof/coverage_exact.json'))
opens=set(); ok=True; info=[]
for nm,(lo,hi) in {'2-48':(2,48),'49-100':(49,100),'101-200':(101,200)}.items():
    c=cov[nm]; sls=[x for x in ALL if lo<=rs(x)[1]<=hi]
    o=set(c['open']); opens|=o
    good = (c['total']==len(sls)) and (c['closed']==len(sls)-len(o)) and o<=set(sls) and len(o)==len(c['open'])
    ok&=good; info.append(f'{nm}: total {c["total"]}/{len(sls)} closed {c["closed"]} open {len(o)}')
chk('1.selmer_coverage_consistent', ok, '; '.join(info))
for sl in ALL:
    if rs(sl)[1]<=200 and sl not in opens:
        add(sl,'А','selmer2isog','criterion_proof/coverage_exact.json',
            how='точный 2-изогенный Селмер семьи E_{m,n} (12 классов), ранг 0; поштучной записи нет — закрытость = дополнение к списку open')
fs=json.load(open(f'{ROOT}/census_genus1/final_status_200.json'))
chk('1.final_status_still_open_empty', fs['still_open']==[], str(fs['still_open']))
chk('1.final_status_keys_eq_selmer_open', set(fs['sources'])==opens, f'sources {len(fs["sources"])}, open {len(opens)}, diff {sorted(set(fs["sources"])^opens)[:10]}')
P('  методы final_status:', dict(collections.Counter(fs['sources'].values())))
six={r['slope']:r for r in jl(f'{ROOT}/census_six_cells/census.jsonl')}
kol={r['slope']:r for r in jl(f'{ROOT}/census_six_cells/kolyvagin.jsonl')}
g1={r['slope']:r for r in jl(f'{ROOT}/census_genus1/g1census.jsonl')}
stage2=jl(f'{ROOT}/magma_calc/stage2.jsonl'); conf15=jl(f'{ROOT}/magma_calc/confirm15.jsonl')
rest3=jl(f'{ROOT}/magma_calc/rest3.jsonl'); last2=jl(f'{ROOT}/magma_calc/last2.jsonl')
def zlist(out):
    m=re.search(r'(?:Z|ZVALS): \[(.*?)\]', out, re.S)
    return None if not m else [x.strip() for x in m.group(1).split(',') if x.strip()]
g2_bad=[]
def g2_ok(sl, zs, cells):
    r,s=rs(sl); L=Lam(r,s)
    if zs is None: return False, 'нет списка z'
    if not set(cells)<=L: return False, f'клетки {cells} не из Λ'
    nd=[z for z in zs if nondeg(z,r,s)]
    return (not nd), (f'невырожденные {nd}' if nd else f'{len(zs)} точек, все вырождены')
mism=[]
for sl,meth in fs['sources'].items():
    r,s=rs(sl); L=Lam(r,s)
    if meth=='шесть клеток (PARI)':
        rec=six.get(sl); good=rec and rec['closed'] and set(rec['by'][0])<=set(abs(x) for x in L) and rec['by'][1] in '+-'
        if not good: mism.append((sl,meth))
        add(sl,'А','six_cells_ellrank','census_six_cells/census.jsonl', how=f'шесть клеток {rec["by"][0]} E{rec["by"][1]}: PARI ellrank верх. 0 + кручение')
    elif meth=='шесть клеток (L(E,1))':
        rec=kol.get(sl); good=rec and rec['closed'] and set(rec['by']['abc'])<=set(abs(x) for x in L)
        if not good: mism.append((sl,meth))
        add(sl,'Б','six_cells_L','census_six_cells/kolyvagin.jsonl', how=f'шесть клеток {rec["by"]["abc"]} E{rec["by"]["kind"]}: L(E,1)/Ω={rec["by"]["L_ratio"]} ≠ 0 (Колывагин)')
    elif meth.startswith('род 1'):
        rec=g1.get(sl); tag='ellrank' if 'ellrank' in meth else 'L_ratio'
        good=rec and rec['closed'] and rec['by'][1]==tag and set(rec['by'][0])<=L
        if not good: mism.append((sl,meth))
        if tag=='ellrank': add(sl,'А','g1_ellrank','census_genus1/g1census.jsonl', how=f'род 1, клетки {rec["by"][0]}: ellrank верх. 0 + кручение')
        else: add(sl,'Б','g1_L','census_genus1/g1census.jsonl', how=f'род 1, клетки {rec["by"][0]}: L(E,1) ≠ 0')
    elif meth=='род 2 (этап 2, Magma)':
        st=[x for x in stage2 if x['slope']==sl and x.get('closed')]
        cf=[x for x in conf15 if x['slope']==sl and 'PROVEN: true' in x['out']]
        good=bool(st) and bool(cf)
        for x in st:
            o,msg=g2_ok(sl,x['z'],[x['c']]+x['pair']); good&=o
        for x in cf:
            o,msg=g2_ok(sl,zlist(x['out']),[int(x['c'])]+[int(y) for y in x['pair']]); good&=o
        if not good: mism.append((sl,meth))
        c=cf[0] if cf else None
        add(sl,'В','g2_magma','magma_calc/confirm15.jsonl; magma_calc/stage2.jsonl',
            how=f'род 2 (5 клеток c={c["c"]}, пара {c["pair"]}): RationalPointsGenus2 флаг true (confirm15) + RankBound=1 и Шаботи (stage2)' if c else 'нет подтверждения')
    elif meth=='род 2 (Claude, Magma)':
        if sl in ('115/163','67/168'):
            x=[x for x in last2 if x['slope']==sl and x.get('closed')]
            good=bool(x) and 'PROVEN: true' in x[0]['out']
            if good:
                o,msg=g2_ok(sl,zlist(x[0]['out']),x[0]['singles']+x[0]['pairs']); good&=o
            if not good: mism.append((sl,meth))
            add(sl,'В','g2_magma','magma_calc/last2.jsonl', how=f'род 2 {x[0]["curve"]}: RationalPointsGenus2 флаг true')
        else:
            xs=[x for x in rest3 if x['slope']==sl and 'PROVEN: true' in x['out']]
            good=bool(xs); msgs=[]
            for x in xs:
                o,msg=g2_ok(sl,zlist(x['out']),[x['single']]+x['pair']); good&=o; msgs.append(msg)
            if not good: mism.append((sl,meth))
            add(sl,'В','g2_magma','magma_calc/rest3.jsonl', how=f'род 2, {len(xs)} кривых (5 клеток) с флагом true; все точки вырождены (переигровка)')
    elif meth=='род 2 (Codex, Magma)':
        txt=open(f'{ROOT}/magma_calc/r_3_86.txt').read()
        good= sl=='3/86' and 'PROVEN: true' in txt
        if good:
            o,msg=g2_ok(sl,zlist(txt),[3,86,89]); good&=o
        if not good: mism.append((sl,meth))
        add(sl,'В','g2_magma','RESULT_3_86_CLOSED_2026-09-13_FROM_CODEX.md; magma_calc/r_3_86.txt',
            how='род 2 y²=(1−3z)(1−86²z²)(1−89²z²): RankBounds=[1,1] + RationalPointsGenus2 флаг true (Codex); повтор Claude — флаг true')
    else:
        mism.append((sl,'неизвестный метод '+meth))
chk('1.final_status_records_match', not mism, f'расхождений {len(mism)}: {mism[:10]}')

# ---------- 2. 201 <= s <= 500
cfiles={'201-300':'census_300/c201_300.jsonl','301-400':'census_300/c301_400.jsonl','401-500':'census_300/c401_500.jsonl'}
crec={}; bad=[]
for nm,fn in cfiles.items():
    R=jl(f'{ROOT}/{fn}'); sls=[x['slope'] for x in R]
    exp=[x for x in ALL if rng(x)==nm]
    chk(f'2.{nm}.complete_no_dups', len(sls)==len(set(sls)) and set(sls)==set(exp), f'{len(sls)} строк, {len(set(sls))} различных, ожидалось {len(exp)}, лишних {len(set(sls)-set(exp))}, недостающих {len(set(exp)-set(sls))}')
    al=[x for x in R if 'ALERT' in x]
    chk(f'2.{nm}.no_ALERT', not al, str(al[:3]))
    for x in R: crec[x['slope']]=(x,fn)
gaps=jl(f'{ROOT}/census_300/c_gaps.jsonl')
for x in gaps:
    sl=x['slope']
    if sl in crec:
        same=crec[sl][0].get('by')==x.get('by')
        notes.append(f'c_gaps: {sl} есть и в {crec[sl][1]} (запись {"совпадает" if same else "ОТЛИЧАЕТСЯ"})')
        if not same: problems.append(f'c_gaps {sl}: запись отличается от {crec[sl][1]}')
    else:
        notes.append(f'c_gaps: {sl} (s ≤ 200) — вне census_300; при s ≤ 200 закрыт основным методом')
        if sl in ALLSET: add(sl,'А','g1_ellrank','census_300/c_gaps.jsonl', how=f'род 1, клетки {x["by"][1]}: ellrank верх. 0 (добор 25.09)', secondary=True)
for sl,(x,fn) in crec.items():
    r,s=rs(sl); L=Lam(r,s)
    if x.get('closed'):
        kind,cells=x['by'][0],x['by'][1]
        if kind=='g1':
            if not set(cells)<=L or len(cells) not in (3,4): bad.append(sl)
            add(sl,'А','g1_ellrank',fn, how=f'род 1, клетки {cells}: PARI ellrank верх. 0 + кручение')
        elif kind=='six':
            if not set(cells)<=set(abs(v) for v in L): bad.append(sl)
            add(sl,'А','six_cells_ellrank',fn, how=f'шесть клеток {cells} E{x["by"][2]}: PARI ellrank верх. 0 + кручение')
        else: bad.append(sl)
chk('2.census_by_cells_in_Lambda', not bad, str(bad[:10]))
unclosed=sorted([sl for sl,(x,fn) in crec.items() if not x.get('closed')], key=key)
P('  не закрыто прогоном census_300:', len(unclosed))
lp={}
for nm in ('201_300','301_400','401_500'):
    for x in json.load(open(f'{ROOT}/census_300/lpass_{nm}.json')):
        if x['slope'] in lp: problems.append(f'lpass дубль {x["slope"]}')
        lp[x['slope']]=(x,f'census_300/lpass_{nm}.json')
chk('2.lpass_covers_unclosed_except_crash', set(lp)==set(unclosed)-{'96/407'}, f'lpass {len(lp)}, unclosed {len(unclosed)}, diff {sorted(set(lp)^set(unclosed))}')
badL=[]
for sl,(x,fn) in lp.items():
    if x.get('closed'):
        r,s=rs(sl); L=Lam(r,s); by=x['by']
        cells=by[1]
        okc = set(cells)<=L if by[0]=='g1' else set(cells)<=set(abs(v) for v in L)
        if not okc or by[-1]!='L_ratio': badL.append(sl)
        # кандидат должен быть из списка lo0 исходной записи
        lo0=crec[sl][0].get('lo0',[])
        if [by[0],by[1]]+([by[2]] if by[0]=='six' else []) not in [c for c in lo0]: badL.append(sl+' (не из lo0)')
        add(sl,'Б','L_ratio',fn, how=f'{"род 1" if by[0]=="g1" else "шесть клеток"}, {by[1]}{(" E"+by[2]) if by[0]=="six" else ""}: L(E,1) ≠ 0 (L_ratio, Колывагин)')
chk('2.lpass_records_consistent', not badL, str(badL))
hard=sorted([sl for sl in unclosed if not (sl in lp and lp[sl][0].get('closed'))], key=key)
P('  трудных после L-прохода:', len(hard), hard)
# 96/407
ex=json.load(open(f'{ROOT}/fable_bl/extra_checks.json'))
e96=ex.get('96/407',[])
okk=[e for e in e96 if e.get('lo')==0 and e.get('hi')==0 and e.get('bad')==[]]
zs_ok=all(not nondeg(z,96,407) for e in okk for z in e['z'])
chk('2.96_407_extra_rank0', bool(okk) and zs_ok and crec['96/407'][0].get('crash'), f'{len(okk)} множителей ellrank 0/0, z вырождены: {zs_ok}')
if okk: add('96/407','А','g1_ellrank_extra','fable_bl/extra_checks.json', how=f'род 1, S={okk[0]["S"]} (k=r/s): ellrank lo=hi=0, все z вырождены; census_300 дал crash (cysignals)')
# Magma род 2 для трудных
mg=[]
for fn in ('last300.jsonl','last400.jsonl','last500.jsonl'):
    for x in jl(f'{ROOT}/magma_calc/{fn}'): mg.append((x,'magma_calc/'+fn))
magma_closed=collections.defaultdict(list); alert_replay=[]
for x,fn in mg:
    sl=x['slope']; r,s=rs(sl)
    if 'PROVEN: true' not in x['out']: continue
    cells=x['singles']+x['pairs']
    zs=zlist(x['out'])
    o,msg=g2_ok(sl,zs,cells)
    # проверка формулы кривой против singles/pairs
    exp_curve=('*'.join([f'(1+{c}*z)' if len(x['singles'])==1 else f'(1+{c}*z)' for c in x['singles']]))
    if x.get('closed') and not o: problems.append(f'{fn} {sl}: closed=true, но переигровка: {msg}')
    if 'ALERT' in x:
        alert_replay.append(dict(slope=sl,file=fn,curve=x['curve'],replay_degenerate=o,msg=msg))
    if x.get('closed') and o:
        magma_closed[sl].append(dict(file=fn,curve=x['curve'],msg=msg))
    if 'ALERT' in x and o:
        magma_closed[sl].append(dict(file=fn,curve=x['curve'],msg=msg,flag='ALERT в записи (старый критерий, коммит e24264e); переигровка: вырождено'))
chk('2.magma_alerts_are_false_alarms', all(a['replay_degenerate'] for a in alert_replay), json.dumps(alert_replay,ensure_ascii=False))
for sl,L_ in magma_closed.items():
    clean=[d for d in L_ if 'flag' not in d]
    if clean:
        d=clean[0]; add(sl,'В','g2_magma',d['file'], how=f'род 2 {d["curve"]}: RationalPointsGenus2 флаг true; {d["msg"]}', n_curves_true=len(clean))
    else:
        d=L_[0]; add(sl,'В','g2_magma',d['file'], how=f'род 2 {d["curve"]}: флаг true; {d["msg"]}', doubt=d['flag'])
# rb_retry: только RankBound (не закрытия)
rbfiles=['magma_calc/rb_retry.jsonl','magma_calc/rb_retry2.jsonl']
rbmin={}
for fn in rbfiles:
    for x in jl(f'{ROOT}/{fn}'):
        m=re.search(r'RB: (\d+)',x['out'])
        if m: rbmin[x['slope']]=min(rbmin.get(x['slope'],99),int(m.group(1)))
notes.append(f'rb_retry*.jsonl: только RankBound для {sorted(rbmin)}; минимальная граница ранга якобиана {rbmin} ≥ 2 ⇒ обычный Шаботи неприменим, это не закрытия')
# Демьяненко ранга 1 (Fable §7)
d45=json.load(open(f'{ROOT}/joint_sieve/dem_open45_safety2.json'))
chk('2.dem45_keys_eq_hard_minus_96_407', set(d45)==set(hard)-{'96/407'}, f'dem {len(d45)}, hard {len(hard)}, diff {sorted(set(d45)^(set(hard)-{"96/407"}))}')
badD=[]
for sl,v in d45.items():
    r,s=rs(sl); L=Lam(r,s)
    okstr= isinstance(v,str) and re.search(r'ЗАКРЫТ \(нет невырожденных решений\): True',v) and 'использовано 3' in v
    fn=f'joint_sieve/dem_{r}_{s}.json'; per=None
    if os.path.exists(f'{ROOT}/{fn}'):
        per=json.load(open(f'{ROOT}/{fn}'))
        res=per['results']; used=[t for t in res if t.get('ok')]
        okper= per.get('closed') is True and len(used)>=3 and all(t['rank']==1 and t['sat_index']==1 and t['nondeg_solutions']==[] and set(t['T'])<=L and set(t['T'])!={-x for x in t['T']} for t in used)
        M0=[t['M0'] for t in used]
    else: okper=False; M0=[]
    if not (okstr and okper): badD.append(sl)
    add(sl,'Г','dem_rank1',f'joint_sieve/dem_open45_safety2.json; {fn}', how=f'Демьяненко–Манин ранга 1 (Fable §7), 3 множителя, M0={M0}, B = Sage silverman_height_bound ×2')
chk('2.dem45_records_consistent', not badD, str(badD))
# малые наклоны с dem-файлами (s<=20 и др.) — вторичное свидетельство
extra_dem=[]
for fn in glob.glob(f'{ROOT}/joint_sieve/dem_*_*.json'):
    b=os.path.basename(fn)
    m=re.fullmatch(r'dem_(\d+)_(\d+)\.json',b)
    if not m: continue
    sl=f'{m.group(1)}/{m.group(2)}'
    if sl in d45: continue
    per=json.load(open(fn))
    if per.get('closed') is True and sl in ALLSET:
        add(sl,'Г','dem_rank1','joint_sieve/'+b, how='Демьяненко–Манин ранга 1 (контроль Fable §7)', secondary=True); extra_dem.append(sl)
notes.append(f'joint_sieve/dem_r_s.json вне 45: {len(extra_dem)} наклонов (контроль §7 на s ≤ 20 и др.), все уже закрыты сильнее')
# строгая B (параллельные чанки rigorous_census/chunk0..2)
strict={}
for ch in ('chunk0','chunk1','chunk2'):
    for fn in sorted(glob.glob(f'{ROOT}/rigorous_census/{ch}/rig_*_*.json')):
        try: d=json.load(open(fn))
        except Exception as e: problems.append(f'не читается {fn}: {e}'); continue
        sl=d.get('slope'); cl=d.get('closed_rigorous', d.get('closed_rigorously'))
        al=d.get('ALERT', d.get('alert'))
        strict[sl]=dict(file=fn.replace(ROOT+'/',''),closed=cl,alert=al,note=os.path.exists(f'{ROOT}/rigorous_census/{ch}/NOTE.md'))
for sl,v in strict.items():
    if v['closed'] is True and not v['alert']:
        add(sl,'Г','dem_rank1_strictB',v['file'], how='Демьяненко–Манин ранга 1 со строгой B (сертификат Безу) — параллельная перепроверка 26.09', note_present=v['note'])
chk('2.strictB_covers_45', set(strict)>=set(d45), f'строгая B: {len(strict)} наклонов, закрыто {sum(1 for v in strict.values() if v["closed"] is True)}, не покрыто {sorted(set(d45)-set(strict))}')
# строгий Демьяненко ранга 2
for sl,fn,extra in (('126/451','rigorous_fable/closure_126_451.json','THEOREM_CLOSURE_126_2026-09-25_FROM_CODEX.md; audit_126/NOTE_AUDIT_126_2026-09-25.md'),
                    ('73/362','rigorous_fable/closure_73_362.json','FABLE_RIGOROUS_CLOSURES_2026-09-25.md; audit_fable_rigorous/NOTE_AUDIT_73_265_2026-09-25.md'),
                    ('265/298','rigorous_fable/closure_265_298.json','FABLE_RIGOROUS_CLOSURES_2026-09-25.md; audit_fable_rigorous/NOTE_AUDIT_73_265_2026-09-25.md')):
    d=json.load(open(f'{ROOT}/{fn}'))
    good = d['slope']==sl and d['case_A_solutions']==[] and d['case_B_solutions']==[]
    chk(f'2.dem2_strict_{sl}', good, f'slope {d["slope"]}, A {d["case_A_solutions"]}, B {d["case_B_solutions"]}, H0 {d.get("H0")}')
    add(sl,'Д','dem_rank2_strict',fn+'; '+extra, how='Демьяненко–Манин ранга 2, константы спаривания доказаны тождествами Безу; 2–3 реализации')
# квадратичный Шаботи — численно
for sl in ('265/298','73/362'):
    add(sl,'Е','qc_bielliptic','FABLE_QC_BIELLIPTIC_2026-09-15.md; qc_bielliptic/fable/', how='квадратичный Шаботи (BP22), «численно с доказательной схемой»', secondary=True)

# ---------- 3. выбор основного основания
# Порядок строгости (обоснование в NOTE): А > Б > В > Д > Г(строгая B) > Г > Е; записи с сомнением (doubt) и вторичные не выбираются основными, если есть чистые.
RANK={('А',None):0,('Б',None):1,('В',None):2,('Д',None):3,('Г','dem_rank1_strictB'):4,('Г',None):5,('Е',None):6}
def rk(e): return RANK.get((e['m'],e['b']), RANK.get((e['m'],None)))
master={}; miss=[]
for sl in ALL:
    E=ev.get(sl,[])
    if not E: miss.append(sl); continue
    prim=[e for e in E if not e.get('doubt') and not e.get('secondary')] or E
    best=min(prim,key=rk)
    alt=[e for e in E if e is not best]
    rec=dict(m=best['m'],b=best['b'],f=best['f'],how=best['how'])
    if best.get('doubt'): rec['doubt']=best['doubt']
    if alt: rec['alt']=[dict((k,v) for k,v in a.items() if k in ('m','b','f','doubt','note_present')) for a in sorted(alt,key=rk)]
    master[sl]=rec
chk('3.every_slope_has_basis', not miss, f'без основания: {len(miss)} {miss[:20]}')
# сводка
tab=collections.OrderedDict()
for nm,a,b in RANGES:
    c=collections.Counter()
    for sl in ALL:
        if rng(sl)==nm and sl in master: c[(master[sl]['m'],master[sl]['b'])]+=1
    tab[nm]={f'{k[0]}:{k[1]}':v for k,v in sorted(c.items())}
P(json.dumps(tab,ensure_ascii=False,indent=1))
weak=[sl for sl in ALL if sl in master and master[sl]['m'] in ('Г','Е')]
weak_detail={sl:dict(primary=master[sl]['b'],alts=[a['m']+':'+a['b'] for a in master[sl].get('alt',[])]) for sl in weak}
meta=dict(
  generated='2026-09-26, Claude (rigorous_census/master/build_master.py); сводка записей, кривые не пересчитывались',
  slopes='r/s, 0<r<s, gcd(r,s)=1, 3<=s<=500 (r≠s/2 автоматически при s≥3); всего 76114',
  legend={'А':'ранг 0: PARI ellrank верхняя граница 0 или точный 2-изогенный Селмер (+ кручение)',
          'Б':'ранг 0: L(E,1)≠0 (L_ratio Sage) + Гросс–Загир–Колывагин, модулярность',
          'В':'род 2: Magma RationalPointsGenus2 с флагом true (онлайн-калькулятор, 14.09)',
          'Г':'Демьяненко–Манин ранга 1 через z↦−z (Fable §7); b=dem_rank1_strictB — перепроверка со строгой B (26.09)',
          'Д':'строгий Демьяненко–Манин ранга 2 (Codex/Fable 25.09, константы спаривания доказаны)',
          'Е':'только численно/условно'},
  strength_order='А > Б > В > Д > Г(strictB) > Г > Е (обоснование — rigorous_census/master/NOTE.md)',
  fields={'m':'класс основания','b':'код метода','f':'файл(ы)-основание','how':'что именно','alt':'другие основания того же наклона','doubt':'сомнение в записи'},
  table=tab, weak=weak_detail)
json.dump(dict(meta=meta,slopes=master),open(OUT_JSON,'w'),ensure_ascii=False,separators=(',',':'))
json.dump(dict(checks=checks,problems=problems,notes=notes,alert_replay=alert_replay,hard=hard,rbmin=rbmin,strict={k:v for k,v in strict.items()},table=tab,weak=weak_detail),open(OUT_CHK,'w'),ensure_ascii=False,indent=1)
P('проблем:',len(problems)); [P('  ',x) for x in problems]
P('заметки:'); [P('  ',x) for x in notes]
P('слабые (основное Г/Е):',len(weak), weak_detail)
P('записано',OUT_JSON,os.path.getsize(OUT_JSON),'байт; время',round(time.time()-t0,1),'с')
