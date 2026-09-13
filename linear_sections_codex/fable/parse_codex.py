"""Собрать из снимка Codex таблицу: k -> {класс набора: границы ранга}, свидетель.
Классы: тройка/четвёрка коэффициентов c, записанных как символы из {1,k,1+k,1-k,-1,-k,-(1+k),-(1-k)},
по модулю c -> -c (замена p -> -p даёт изоморфную кривую)."""
import json,sys
from fractions import Fraction as F
from collections import Counter,defaultdict
D='/home/kep/magicKube/linear_sections_codex/'
def sym(c,k):
    m={1:'1',k:'k',1+k:'1+k',1-k:'1-k',-1:'-1',-k:'-k',-1-k:'-(1+k)',k-1:'-(1-k)'}
    return m[c]
NEG={'1':'-1','k':'-k','1+k':'-(1+k)','1-k':'-(1-k)'}
NEG.update({v:u for u,v in NEG.items()})
ORDER=['1','k','1+k','1-k','-1','-k','-(1+k)','-(1-k)']
def canon(syms):
    a=tuple(sorted(syms,key=ORDER.index)); b=tuple(sorted((NEG[s] for s in syms),key=ORDER.index))
    return min(a,b,key=lambda t:[ORDER.index(x) for x in t])
census=json.load(open(D+'census_through_48.json'))
files={}
for fn in ['extended_slopes_7_12.json','extended_slopes_13_24.json','extended_slopes_25_48.json','four_cell_25_48.json','four_cell_followup.json']:
    files[fn]=json.load(open(D+fn))
small=json.load(open(D+'small_slopes_verified.json'))
table={}   # k -> dict(kind, cls) -> rank_bounds (list or None)
witness={} # k -> (method, cls, curve_key)
for fn,d in files.items():
    rc=d['rank_cache']
    for row in d['slopes']:
        k=F(row['k']); t=table.setdefault(k,{})
        for i,tr in enumerate(row['tried']):
            cs=[F(c) for c in tr['coefficients']]
            cls=canon([sym(c,k) for c in cs])
            kind='C' if len(cs)==3 else 'Q'
            rb=rc[tr['curve_key']].get('rank_bounds')
            key=(kind,cls)
            if key in t and t[key]!=rb and rb is not None and t[key] is not None:
                print('conflict',k,key,t[key],rb,file=sys.stderr)
            if key not in t or t[key] is None: t[key]=rb
            if row.get('status')=='excluded' and row.get('witness_index')==i:
                witness[k]=(kind,cls,tr['curve_key'])
for w in small['rank_zero_exclusions']:
    k=F(w['k']); cs=[F(c) for c in w['selected_coefficients']]
    witness[k]=('C',canon([sym(c,k) for c in cs]),','.join(w['minimal_ainvs']))
    table.setdefault(k,{})[('C',witness[k][1])]=[0,0]
exc=json.load(open(D+'exception_40_43_certified.json'))
out={'table':{str(k):{kind+'|'+','.join(cls):rb for (kind,cls),rb in t.items()} for k,t in table.items()},
     'witness':{str(k):{'kind':v[0],'cls':','.join(v[1]),'curve_key':v[2]} for k,v in witness.items()}}
json.dump(out,open(D+'fable/codex_table.json','w'),indent=0)
print('slopes in table',len(table),'witnesses',len(witness))
print('census methods',Counter(w['method'] for w in census['witnesses']))
# consistency with census
for w in census['witnesses']:
    k=F(w['k'])
    if w['method'] in('cubic','quartic'):
        assert k in witness, k
        if 'curve_key' in w: assert witness[k][2]==w['curve_key'],(k,witness[k],w)
k=F('40/43'); witness[k]=('Q',canon(['k','1+k','-k','-(1+k)']),'0,1,0,-68110640,204180681300'); table[k][('Q',witness[k][1])]=[0,0]
out={'table':{str(k):{kind+'|'+','.join(cls):rb for (kind,cls),rb in t.items()} for k,t in table.items()},
     'witness':{str(k):{'kind':v[0],'cls':','.join(v[1]),'curve_key':v[2]} for k,v in witness.items()}}
json.dump(out,open(D+'fable/codex_table.json','w'),indent=0)
print('census curve_keys match witness table')
