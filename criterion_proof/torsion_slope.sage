# Для наклона k: кривая Y^2 = (1+p)(1+kp)(1+(1+k)p) ~ E_{1,k}. Перечисляем 8 точек кручения ℤ/2×ℤ/4
# как сечения над Q(k), переводим в значения p и смотрим, что происходит с восемью клетками 1 + c p.
K.<k> = FunctionField(QQ)
A = k*(1+k); B = 1 + k + (1+k); C = 1 + k*1 + (1+k) + k*(1+k)   # коэффициенты f(p)=A p^3 + B' p^2 + C' p + 1
R.<p> = K[]
f = (1+p)*(1+k*p)*(1+(1+k)*p)
A,Bc,Cc = f[3],f[2],f[1]
E = EllipticCurve(K,[0,Bc,0,A*Cc,A^2])          # X = A p, V = A Y
cells = {'1+p':1+p,'1-(1+k)p':1-(1+k)*p,'1+kp':1+k*p,'1-(1-k)p':1-(1-k)*p,
         '1+(1-k)p':1+(1-k)*p,'1-kp':1-k*p,'1+(1+k)p':1+(1+k)*p,'1-p':1-p}
P0 = E(0, A)                                     # p = 0
T = [E(r,0) for r in (-(A/1), -(A/k), -(A/(1+k)))]   # корни X = -A/c
pts = [E(0), P0, -P0] + T + [P0+t for t in T]
seen=set()
for P in pts:
    if P in seen: continue
    seen.add(P); seen.add(-P) if P!=E(0) else None
    if P == E(0): print("O  → p = ∞ (вырожденная точка на бесконечности)"); continue
    X,V = P.xy(); pv = X/A
    vals = {nm: c(pv) for nm,c in cells.items()}
    print(f"точка порядка {P.order(algorithm='generic_small')}: p = {pv}")
    for nm,v in vals.items(): print(f"      {nm:10} = {v}")
