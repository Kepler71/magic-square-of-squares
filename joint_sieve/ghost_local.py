from jsieve import *
f = Factor([-18,-17,-1]); print(f.name, 'ранг', f.rank, 'кручение', f.tinv, 'P0', f.decompose(f.point(0,1)))
G = f.gens[0]; t1 = f.tgens[0]; assert t1.order() == 4
D1 = 24*G; D2 = 24*G - 2*t1
res = {'24G∈48E': 0, '24G-2t1∈48E': 0, 'ни то ни другое': 0, 'оба': 0}; bad = []
for l in prime_range(5, 2000):
    if l in f.bad: continue
    El = f.E.change_ring(GF(l)); Gr = El.abelian_group()
    inv = [int(g.order()) for g in Gr.gens()]
    def in48(P):
        v = [int(x) for x in Gr.discrete_log(El(P))] if not P.is_zero() else [0]*len(inv)
        return all(x % gcd(48, m) == 0 for x, m in zip(v, inv))
    a, b = in48(D1), in48(D2)
    if a and b: res['оба'] += 1
    elif a: res['24G∈48E'] += 1
    elif b: res['24G-2t1∈48E'] += 1
    else: res['ни то ни другое'] += 1; bad.append(int(l))
print(res, 'нарушения:', bad[:10])
# глобально: 24G ∈ 48E(ℚ)? 24G−2t1 ∈ 48E(ℚ)? (E(ℚ)=ℤG⊕T ⇒ 48E(ℚ)=48ℤG) — нет. Делимость на 2,4,8 в E(ℚ):
for k in (2,4,8,16):
    print(f'24G делится на {k} в E(ℚ):', len((24*G).division_points(k))>0, f'; 24G−2t1 делится на {k}:', len((24*G-2*t1).division_points(k))>0)
