# Структурная оценка: для каждого набора клеток S (|S| = 3 или 4) на прямой q = k p
# выясняем, является ли вырожденная точка p = 0 точкой кручения на кривой Y^2 = prod(1 + c_i p).
# Если точка бесконечного порядка при случайных k, то (для неизотривиального семейства) по Сильвермену
# ранг >= 1 почти для всех k, и набор S может закрыть лишь конечное число наклонов.
import itertools, json
R.<kk> = QQ[]
names = ["1+p","1-(1+k)p","1+kp","1-(1-k)p","1+(1-k)p","1-kp","1+(1+k)p","1-p"]
coef  = [R(1), -(1+kk), kk, -(1-kk), (1-kk), -kk, (1+kk), R(-1)]

def curve_and_point(cs, kv):
    cs = [QQ(c(kv)) for c in cs]
    if len(set(cs)) < len(cs) or 0 in cs: return None, None
    Pp.<p> = QQ[]
    f = prod(1 + c*p for c in cs)
    if len(cs) == 3:
        A, B, C = f[3], f[2], f[1]
        E = EllipticCurve(QQ, [0, B, 0, A*C, A^2])
        P = E(0, A)                       # p = 0
    else:
        c1 = cs[0]; r = -1/c1
        # x = 1/(p - r), w = Y/(p - r)^2 ; w^2 = c1 * prod_{i>=2} ((1 - c_i/c1) x + c_i)
        Px.<x> = QQ[]
        g = c1 * prod((1 - c/c1)*x + c for c in cs[1:])
        A, B, C, D = g[3], g[2], g[1], g[0]
        E = EllipticCurve(QQ, [0, B, 0, A*C, A^2*D])
        x0 = 1/(0 - r); w0 = 1/(0 - r)^2
        P = E(A*x0, A*w0)
    return E, P

ks = [QQ(3)/7, QQ(5)/11, QQ(7)/13, QQ(2)/9, QQ(11)/17, QQ(13)/29]
out = []
for size in (3, 4):
    for S in itertools.combinations(range(8), size):
        orders, js = [], set()
        for kv in ks:
            E, P = curve_and_point([coef[i] for i in S], kv)
            if E is None: orders.append(None); continue
            o = P.order()
            orders.append(str(o))
            js.add(E.j_invariant())
        tors = all(o is not None and o != '+Infinity' for o in orders)
        out.append(dict(S=[names[i] for i in S], idx=list(S), orders=orders, torsion_all=tors, j_const=(len(js)==1)))
json.dump(out, open("p0_torsion.json","w"), ensure_ascii=False, indent=1)
for size in (3,4):
    rows = [o for o in out if len(o['idx'])==size]
    t = [o for o in rows if o['torsion_all']]
    print(f"|S|={size}: наборов {len(rows)}, p=0 кручение при всех k: {len(t)}, j постоянен: {sum(o['j_const'] for o in rows)}")
    for o in t: print("   кручение:", o['idx'], o['S'], "порядки", o['orders'])
