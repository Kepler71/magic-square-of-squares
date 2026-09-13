import json, itertools, random
load("gen_sections.sage") if False else None
R.<kk> = QQ[]
coef = [R(1), -(1+kk), kk, -(1-kk), (1-kk), -kk, (1+kk), R(-1)]
data = json.load(open("p0_torsion.json"))
torsion_sets = {tuple(o['idx']) for o in data if o['torsion_all']}
inf_seen = {tuple(o['idx']): o['orders'].count('+Infinity') for o in data if not o['torsion_all']}
print("наборов без кручения, у которых бесконечный порядок встретился хотя бы раз:", sum(1 for v in inf_seen.values() if v>0), "из", len(inf_seen))

def cubic_curve(S, kv):
    cs = [QQ(coef[i](kv)) for i in S]
    Pp.<p> = QQ[]; f = prod(1 + c*p for c in cs)
    A,B,C = f[3],f[2],f[1]
    return EllipticCurve(QQ,[0,B,0,A*C,A^2])

cen = json.load(open("/home/kep/magicKube/linear_sections_codex/census_through_48.json"))
cub = [w for w in cen['witnesses'] if w['method']=='cubic']
random.seed(int(20260913)); sample = random.sample(cub, int(60))
hit_t = hit_n = nomatch = 0; notes=[]
for w in sample:
    kv = QQ(w['k']); Ew = EllipticCurve(QQ,[QQ(a) for a in w['curve_key'].split(',')])
    match = [S for S in itertools.combinations(range(8),3)
             if len(set(QQ(coef[i](kv)) for i in S))==3 and cubic_curve(S,kv).is_isomorphic(Ew)]
    if not match:
        # допускаем изогению/твист модели: сравним по j
        match = [S for S in itertools.combinations(range(8),3)
                 if len(set(QQ(coef[i](kv)) for i in S))==3 and cubic_curve(S,kv).j_invariant()==Ew.j_invariant()]
        if match: notes.append((w['k'],'только по j'))
    if not match: nomatch += 1; continue
    if any(tuple(S) in torsion_sets for S in match): hit_t += 1
    else: hit_n += 1; notes.append((w['k'], [list(S) for S in match]))
print(f"выборка кубических свидетелей: {len(sample)}; из 12 «кручёных» троек: {hit_t}; из остальных: {hit_n}; не сопоставлено: {nomatch}")
for n in notes[:10]: print("   ", n)

# знаки функционального уравнения у 22 «хороших» кривых на случайных k
good = [o for o in data if o['torsion_all']]
ks = sorted({QQ(a)/b for b in range(50,90) for a in range(1,b) if gcd(a,b)==1})
random.seed(int(7)); ks = random.sample(ks, int(40))
print("\nзнак ф.у. (+1 долей) по 40 случайным k со знаменателем 50..89:")
for o in good:
    S = o['idx']; signs=[]
    for kv in ks:
        cs=[QQ(coef[i](kv)) for i in S]
        if len(set(cs))<len(cs): continue
        if len(S)==3: E = cubic_curve(S,kv)
        else:
            Px.<x>=QQ[]; c1=cs[0]; g=c1*prod((1-c/c1)*x+c for c in cs[1:])
            A,B,C,D=g[3],g[2],g[1],g[0]; E=EllipticCurve(QQ,[0,B,0,A*C,A^2*D])
        signs.append(E.root_number())
    print(f"   {str(S):16} +1: {signs.count(1):2}/{len(signs)}   кручение при k={ks[0]}: {E.torsion_subgroup().invariants()}")
