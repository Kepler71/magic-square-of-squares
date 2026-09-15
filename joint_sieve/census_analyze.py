# Разбор переписи совместного решета
import json, sys, collections
fn = sys.argv[1] if len(sys.argv) > 1 else 'census_3_40.jsonl'
rows = [json.loads(l) for l in open(fn)]
ok = [r for r in rows if 'err' not in r]
err = [r for r in rows if 'err' in r]
print('всего', len(rows), 'посчитано', len(ok), 'ошибок/пропусков', len(err))
print('пропуски:', collections.Counter(r['err'] for r in err))
ex = [r for r in ok if r['extra'] > 0]
miss = [r for r in ok if r['missing'] > 0]
print('с невырожденными выжившими:', len(ex), [(r['slope'], r['extra']) for r in ex][:20])
print('с потерянными известными (ошибка!):', len(miss))
print('вырожденных z помимо 0:', [(r['slope'], r['degz']) for r in ok if len(r['degz']) > 1][:10])
print('выживших = известных:', sum(1 for r in ok if r['survivors'] == r['known']), 'из', len(ok))
print('распределение известных классов:', collections.Counter(r['known'] for r in ok))
print('наклоны с множителем ранга 0 (закрыты рангом 0):', sum(1 for r in ok if r['rank0'] > 0), '; без ранга 0:', sum(1 for r in ok if r['rank0'] == 0))
print('без ранга 0 и с невырожденными:', [(r['slope'], r['extra']) for r in ok if r['rank0'] == 0 and r['extra'] > 0])
# P_T(0): образующая?
gen = collections.Counter()
for r in ok:
    for rk, c, t in zip(r['ranks'], r['p0'], r['tors']):
        n = c[:rk]
        if all(x == 0 for x in n): gen['кручение'] += 1
        elif rk == 1 and abs(n[0]) == 1: gen['±образующая (r=1)'] += 1
        elif rk == 1: gen[f'кратная (r=1), n={abs(n[0])}'] += 1
        else: gen['r=2: примитивная' if any(abs(x) == 1 for x in n) else 'r=2: другое'] += 1
print('P_T(0) в использованных множителях:', dict(gen))
print('индексы насыщения >1:', collections.Counter(i for r in ok for i in r['sat_index'] if i > 1))
print('макс. кандидатов:', max(r['maxcand'] for r in ok), 'среднее время', sum(r['time'] for r in ok) / len(ok))
print('простых в среднем', sum(r['nprimes'] for r in ok) / len(ok))
