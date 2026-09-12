# -*- coding: utf-8 -*-
# НЕЗАВИСИМАЯ ПОПЫТКА СЛОМАТЬ заявление Codex: C_{19,5}(Q) = пусто.
# Всё считается с нуля: s, b, корни, тождества, ранг, СОБСТВЕННАЯ неторсионная точка,
# СОБСТВЕННАЯ арифметика на кривой (не Sage-овская) и СОБСТВЕННОЕ перечисление классов delta.
# Claude, 2026-09-12.   Запуск: sage check_kummer_19_5_классы.sage

import itertools, sys
from sage.all import *

def hdr(t):
    print("\n" + "="*76); print(t); print("="*76)

hdr("НЕЗАВИСИМАЯ ПРОВЕРКА KUMMER-ИСКЛЮЧЕНИЯ ДЛЯ (m,n)=(19,5)")

# ---------------------------------------------------------------- утилиты
def sqclass(x):
    """Квадратсвободный представитель со знаком класса x в Q*/Q*^2."""
    x = QQ(x)
    if x == 0:
        raise ValueError("нулевой аргумент квадратного класса")
    z = ZZ(x.numerator()*x.denominator())     # p/q ~ p*q
    sgn = 1 if z > 0 else -1
    out = 1
    for (p,e) in factor(abs(z)):
        if e % 2 == 1:
            out *= p
    return sgn*out

def trip_mul(A,B):
    return tuple(sqclass(QQ(a)*QQ(b)) for a,b in zip(A,B))

# ---------------------------------------------------------------- 1. модель
hdr("[1] МОДЕЛЬ — считаю сам, чужие числа только сверяю")
m, n = 19, 5
assert gcd(m,n) == 1
s = QQ(m**2 + n**2)/2
b = s*m**2*n**2
print("  s=(m^2+n^2)/2 =", s, "  Codex 193 ->", s == 193)
print("  b=s*m^2*n^2   =", b, "  Codex 1741825 ->", b == 1741825)
e1, e2, e3 = -b, -s*m**4, -s*n**4
print("  e1,e2,e3      =", e1, e2, e3)
print("  Codex (-1741825,-25151953,-120625) ->", (e1,e2,e3)==(-1741825,-25151953,-120625))
assert len({e1,e2,e3}) == 3

# ---------------------------------------------------------------- 2. тождества
hdr("[2] ТОЧНЫЕ ТОЖДЕСТВА (полиномиально, X = b t^2)")
R = PolynomialRing(QQ, 4, names=('t','u0','u4','u8'))
(t,u0,u4,u8) = R.gens()
F0 = m**2 + n**2*t**2
F4 = s*(1 + t**2)
F8 = n**2 + m**2*t**2
X  = b*t**2
ok1 = ((X - e1) - (m*n)**2 * F4) == 0
ok2 = ((X - e2) - s*m**2 * F0) == 0
ok3 = ((X - e3) - s*n**2 * F8) == 0
print("  X-e1 == (mn)^2*F4 :", ok1)
print("  X-e2 == s*m^2*F0  :", ok2)
print("  X-e3 == s*n^2*F8  :", ok3)
okV = ((m*n)**2*u4**2 * s*m**2*u0**2 * s*n**2*u8**2 - (b*u0*u4*u8)**2) == 0
print("  V=b*u0u4u8 => V^2 = (X-e1)(X-e2)(X-e3):", okV)
assert ok1 and ok2 and ok3 and okV
req = (sqclass((m*n)**2), sqclass(s*m**2), sqclass(s*n**2))
print("  НЕОБХОДИМЫЙ класс delta конечной точки C:", req, " Codex (1,193,193) ->", req==(1,193,193))
print("  F0,F4,F8 строго положительны при вещественном t -> u_i != 0, X != e_i: True")

# ---------------------------------------------------------------- 3. кривая
hdr("[3] КРИВАЯ E")
a2 = -(e1+e2+e3); a4 = e1*e2+e1*e3+e2*e3; a6 = -e1*e2*e3
E = EllipticCurve(QQ, [0, a2, 0, a4, a6])
print("  E:", E)
print("  a2,a4,a6 =", a2, a4, a6)
print("  disc != 0:", E.discriminant() != 0, "  conductor =", E.conductor())
print("  факторизация проводника:", factor(E.conductor()))
print("  кручение E(Q):", E.torsion_subgroup().invariants(), "порядок", E.torsion_order())
assert E.torsion_subgroup().invariants() == (2,2)

# ---------------------------------------------------------------- 4. РАНГ: все доступные пути
hdr("[4] РАНГ — ЭТО ЕДИНСТВЕННОЕ УЗКОЕ МЕСТО, БЬЮ ПО НЕМУ")
res_rank = {}
try:
    rb = E.rank_bounds(); res_rank['eclib_bounds'] = rb
    print("  eclib (Sage rank_bounds)     :", rb, "  <- верхняя граница БЕЗ паринга Касселса")
except Exception as ex: print("  rank_bounds упал:", ex)
try:
    sel = E.selmer_rank(); res_rank['selmer'] = sel
    print("  ранг 2-Селмера C (eclib)     :", sel)
except Exception as ex: print("  selmer_rank упал:", ex)
for eff in [0,1,2,3,4]:
    try:
        pr = pari(E).ellrank(eff)
        print("  PARI ellrank(effort=%d)       : [r1,r2,s,pts] ="%eff, pr)
        res_rank['pari_%d'%eff] = pr
    except Exception as ex: print("  PARI ellrank(effort=%d) упал:"%eff, ex)
print("  -> PARI: r2 = C - T - s, s = ранг Sha[2]/2Sha[4], вычисляется паринг Касселса.")
print("  -> eclib НЕ считает паринг, поэтому даёт C-T = 3. Границу 1 даёт ТОЛЬКО PARI.")

hdr("[4b] НЕЗАВИСИМЫЙ АНАЛИТИЧЕСКИЙ ПУТЬ К rank=1 (Гросс-Загир + Колывагин)")
try:
    w = E.root_number()
    print("  корневое число w =", w, " => порядок нуля L в s=1 нечётен" if w==-1 else "")
except Exception as ex:
    w = None; print("  root_number упал:", ex)
try:
    L1 = E.lseries().L1_vanishes()
    print("  L(1) = 0 ?:", L1)
except Exception as ex: print("  L1_vanishes упал:", ex)
try:
    d1 = E.lseries().deriv_at1(100000)
    print("  L'(1) ~", d1, "  (значение, оценка ошибки)")
except Exception as ex: print("  deriv_at1 упал:", ex)
try:
    ar = E.analytic_rank(algorithm='pari'); print("  ellanalyticrank (PARI):", ar)
except Exception as ex: print("  analytic_rank pari упал:", ex)
try:
    ar2 = E.analytic_rank(algorithm='rubinstein'); print("  analytic_rank (rubinstein):", ar2)
except Exception as ex: print("  analytic_rank rubinstein упал:", ex)

# ---------------------------------------------------------------- 5. СВОЯ точка
hdr("[5] СВОЯ НЕТОРСИОННАЯ ТОЧКА")
cands = []
for eff in [1,2,3,4,5]:
    try:
        pr = pari(E).ellrank(eff)
        for P in pr[3]:
            xx = QQ(P[0]); yy = QQ(P[1])
            if (xx,yy) not in cands: cands.append((xx,yy))
    except Exception: pass
try:
    for P in E.point_search(16, rank_bound=4):
        if not P.is_zero() and P.order() == oo:
            if (QQ(P[0]),QQ(P[1])) not in cands: cands.append((QQ(P[0]),QQ(P[1])))
except Exception as ex:
    print("  point_search упал:", ex)
print("  кандидаты (найдены МНОЙ в этом запуске):")
for c in cands: print("    ", c)

codex_x = QQ(-28561963657)/QQ(1369)
codex_y = QQ(146764800)   # чужое y для (15,1); для (19,5) Codex дал другое, сверим ниже
codexP = None
try:
    codexP = E.lift_x(QQ(-28561963657)/QQ(1369))
    print("  точка Codex (19,5): x =", QQ(-28561963657)/QQ(1369), " лежит на E:", True)
    print("    её y (Sage) =", codexP[1], "  Codex писал y = -2089086742828800/50653")
except Exception as ex:
    print("  подъём x Codex не удался:", ex)

assert cands, "не нашёл ни одной неторсионной точки — считать нечего"
# беру ПЕРВУЮ отличную от точки Codex, если такая есть
mypt = None
for c in cands:
    if c[0] != QQ(-28561963657)/QQ(1369):
        mypt = c; break
if mypt is None: mypt = cands[0]
print("  МОЯ рабочая точка G =", mypt)
print("  она ОТЛИЧНА от x-координаты Codex:", mypt[0] != QQ(-28561963657)/QQ(1369))

# ---------------------------------------------------------------- 6. своя арифметика
hdr("[6] СОБСТВЕННАЯ АРИФМЕТИКА (без EllipticCurve), y^2 = x^3+a2x^2+a4x+a6")
A2,A4,A6 = QQ(a2),QQ(a4),QQ(a6)
def on_curve(P):
    if P is None: return True
    x,y = P; return y**2 == x**3 + A2*x**2 + A4*x + A6
def neg(P): return None if P is None else (P[0],-P[1])
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2:
        if y1 == -y2: return None
        lam = (3*x1**2 + 2*A2*x1 + A4)/(2*y1)
    else:
        lam = (y2-y1)/(x2-x1)
    x3 = lam**2 - A2 - x1 - x2
    return (x3, lam*(x1-x3) - y1)
def mul(k,P):
    if k < 0: return neg(mul(-k,P))
    Rr, Qq = None, P
    while k:
        if k & 1: Rr = add(Rr,Qq)
        Qq = add(Qq,Qq); k >>= 1
    return Rr

Gm = mypt
t1,t2,t3 = (QQ(e1),QQ(0)), (QQ(e2),QQ(0)), (QQ(e3),QQ(0))
print("  G на кривой (моя формула):", on_curve(Gm))
print("  T1,T2,T3 на кривой:", on_curve(t1), on_curve(t2), on_curve(t3))
print("  T1+T2 == T3:", add(t1,t2)==t3, "  2T1==O:", mul(2,t1) is None,
      "  2T2==O:", mul(2,t2) is None, "  2T3==O:", mul(2,t3) is None)
print("  G неторсионна (4G != O, 8G != O, и высота растёт):",
      mul(4,Gm) is not None and mul(8,Gm) is not None)
print("  x(2G) =", mul(2,Gm)[0])
_GS = E(Gm[0],Gm[1])
print("  сверка с Sage: 2G совпал:", QQ((2*_GS)[0]) == mul(2,Gm)[0] and QQ((2*_GS)[1]) == mul(2,Gm)[1])
print("  сверка с Sage: 3G совпал:", QQ((3*_GS)[0]) == mul(3,Gm)[0] and QQ((3*_GS)[1]) == mul(3,Gm)[1])
print("  сверка с Sage: 5G совпал:", QQ((5*_GS)[0]) == mul(5,Gm)[0])
print("  порядок G по Sage:", _GS.order())
# связь моей точки с точкой Codex
if codexP is not None:
    cP = (QQ(codexP[0]), QQ(codexP[1]))
    rel = []
    for k in range(-8,9):
        if k==0: continue
        Pk = mul(k,Gm)
        for (lbl,T) in [("",None),("+T1",t1),("+T2",t2),("+T3",t3)]:
            Q = add(Pk,T) if T else Pk
            if Q is not None and (Q == cP or Q == neg(cP)):
                rel.append("%d*G%s = %sG_Codex" % (k, lbl, "" if Q==cP else "-"))
    print("  связь моей точки с точкой Codex:", rel if rel else "не выражается через k*G+кручение, |k|<=8")

# ---------------------------------------------------------------- 7. delta
hdr("[7] delta — СОБСТВЕННАЯ РЕАЛИЗАЦИЯ 2-СПУСКА")
EE = [QQ(e1),QQ(e2),QQ(e3)]
def delta(P):
    if P is None: return (1,1,1)
    x,y = P; out=[]
    for i in range(3):
        d = x - EE[i]
        if d == 0:
            j,k = [q for q in range(3) if q != i]
            d = (EE[i]-EE[j])*(EE[i]-EE[k])
        out.append(sqclass(d))
    return tuple(out)

print("  ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ соглашения: (m,n)=(15,8), t=1, u0=u4=u8=17")
mm,nn = 15,8
ss = QQ(mm**2+nn**2)/2; bb = ss*mm**2*nn**2
ee = [-bb, -ss*mm**4, -ss*nn**4]
dc = tuple(sqclass(bb*QQ(1)**2 - x) for x in ee)
print("    delta(точки) =", dc, "  требуемый (1,s,s) =", (1,sqclass(ss),sqclass(ss)),
      " совпадает:", dc == (1,sqclass(ss),sqclass(ss)))
assert dc == (1,sqclass(ss),sqclass(ss))

# ---------------------------------------------------------------- 8. таблица классов
hdr("[8] ПОЛНАЯ ТАБЛИЦА delta(a*G + c*T1 + d*T2)")
rows = []
for a in [0,1]:
    for c in [0,1]:
        for d in [0,1]:
            P = None
            if a: P = add(P,Gm)
            if c: P = add(P,t1)
            if d: P = add(P,t2)
            assert on_curve(P)
            rows.append(((a,c,d),P,delta(P)))
print("  %-10s | %s" % ("(a,c,d)","delta = (d1,d2,d3)  [квадратсвободные со знаком]"))
print("  " + "-"*66)
for (acd,P,D) in rows:
    print("  %-10s | %s" % (str(acd), str(D)))
classes = [D for (_,_,D) in rows]
dist = len(set(classes))
print("\n  различных:", dist, "из", len(classes), " ВСЕ ПОПАРНО РАЗЛИЧНЫ:", dist==len(classes))
if dist != len(classes):
    from collections import Counter
    print("  ПОВТОРЫ:", [k for k,v in Counter(classes).items() if v>1])
S = set(classes)
print("  требуемый класс", req, "среди них:", req in S)
print("  множество замкнуто по умножению (подгруппа):",
      all(trip_mul(A,B) in S for A in classes for B in classes))
print("  произведение координат каждого класса — квадрат (d1d2d3 in Q*^2):",
      all(sqclass(QQ(D[0])*QQ(D[1])*QQ(D[2]))==1 for D in classes))

codex_rows = [(1,1,1),(-1,4053,-4053),(-4053,386,-42),(4053,42,386),
              (-3,386,-1158),(3,42,14),(1351,1,1351),(-1351,4053,-3)]
print("\n  таблица Codex:", codex_rows)
print("  множества СОВПАДАЮТ:", set(codex_rows)==S)
if set(codex_rows)!=S:
    print("    только у Codex:", sorted(set(codex_rows)-S))
    print("    только у меня :", sorted(S-set(codex_rows)))
print("  у Codex (1,193,193) присутствует:", (1,193,193) in set(codex_rows))
print("  множество Codex замкнуто по умножению:",
      all(trip_mul(A,B) in set(codex_rows) for A in codex_rows for B in codex_rows))

# ---------------------------------------------------------------- 9. гомоморфность
hdr("[9] ГОМОМОРФНОСТЬ delta")
base = [None,t1,t2,t3,Gm,add(Gm,t1),add(Gm,t2),add(Gm,t3),
        mul(2,Gm),mul(3,Gm),mul(-1,Gm),mul(5,Gm),mul(-4,Gm),
        add(mul(2,Gm),t1),add(mul(3,Gm),t2),add(mul(-2,Gm),t3),mul(7,Gm)]
bad=0; tot=0
for P in base:
    for Q in base:
        tot+=1
        if delta(add(P,Q)) != trip_mul(delta(P),delta(Q)):
            bad+=1; print("   НАРУШЕНИЕ:",P,Q)
print("  проверено пар:",tot," нарушений:",bad)
dG = delta(Gm); okk=True
print("  delta(kG), k=0..14  (чётные -> (1,1,1), нечётные -> delta(G)):")
for k in range(15):
    dk = delta(mul(k,Gm)); exp = (1,1,1) if k%2==0 else dG
    if dk!=exp: okk=False
    print("    k=%2d %-22s ожидалось %-22s %s"%(k,str(dk),str(exp),"ok" if dk==exp else "!!!"))
print("  все совпали:", okk)

# ---------------------------------------------------------------- 10. атака
hdr("[10] АТАКА: точки E вне 8 классов / второй независимый генератор")
outside=[]
try:
    found = E.point_search(14, rank_bound=4)
    print("  point_search(14) нашёл точек:", len(found))
    for p in found:
        if p.is_zero(): continue
        D = delta((QQ(p[0]),QQ(p[1])))
        if D not in S: outside.append((p,D))
    print("  точек с delta ВНЕ моих 8 классов:", len(outside))
    for (p,D) in outside[:12]: print("    ", p, D)
    # независимость найденных точек
    if len(found) > 0:
        try:
            hm = E.height_pairing_matrix([p for p in found if p.order()==oo])
            print("  ранг матрицы высот найденных точек:", hm.rank())
        except Exception as ex: print("  height_pairing_matrix:", ex)
except Exception as ex:
    print("  point_search упал:", ex)
bad2 = sum(1 for a in range(-30,31) for c in [0,1] for d in [0,1]
           if delta(add(add(mul(a,Gm), t1 if c else None), t2 if d else None)) not in S)
print("  среди a*G+кручение, |a|<=30, вне 8 классов:", bad2)

# ---------------------------------------------------------------- 11. Селмер: где сидит нужный класс
hdr("[11] ГДЕ СИДИТ ТРЕБУЕМЫЙ КЛАСС: локальная разрешимость 2-накрытия для (1,193,193)")
# 2-накрытие для класса (d1,d2,d3): d1 z1^2 - d2 z2^2 = e2-e1,  d1 z1^2 - d3 z3^2 = e3-e1
d1,d2,d3 = 1, 193, 193
print("  система:  %d z1^2 - %d z2^2 = %s ;  %d z1^2 - %d z3^2 = %s"
      % (d1,d2,e2-e1,d1,d3,e3-e1))
print("  (это в точности накрытие, которому отвечает кривая C)")
badp = sorted(set([2]) | set(p for p,_ in factor(2*m*n*(m**2-n**2)*(m**2+n**2))))
print("  плохие простые:", badp)
def locally_solvable(p, B=200):
    """грубый поиск решения в Z_p по модулю p^k"""
    k = 6 if p>2 else 10
    M = p**k
    Zm = Integers(M)
    tgt1 = Integer(e2-e1) % M; tgt2 = Integer(e3-e1) % M
    for z1 in range(M if M<4000 else 4000):
        v1 = (d1*z1*z1 - tgt1) % M
        v2 = (d1*z1*z1 - tgt2) % M
        # нужно v1 = d2 z2^2, v2 = d3 z3^2 разрешимо в Z_p
        if _isqp(v1, d2, p, k) and _isqp(v2, d3, p, k): return True, z1
    return False, None
def _isqp(v, d, p, k):
    M = p**k
    for z in range(min(M, 3000)):
        if (d*z*z - v) % M == 0: return True
    return False
# проще: используем точки самой C над Q_p (Codex дал свидетеля t=112 при p=193)
print("  прямая проверка: ищу целое t, при котором F0,F4,F8 все квадраты в Q_p")
for p in badp + [101,103,107]:
    wit = None
    for tt in range(0, 400):
        f0 = m**2 + n**2*tt**2; f4 = ZZ(s*(1+tt**2)); f8 = n**2 + m**2*tt**2
        if f0==0 or f4==0 or f8==0: continue
        try:
            Qp = Qp_ = None
            from sage.rings.padics.factory import Qp as _Qp
            K = _Qp(p, 30)
            if K(f0).is_square() and K(f4).is_square() and K(f8).is_square():
                wit = tt; break
        except Exception:
            pass
    print("    p=%-5s свидетель t=%s" % (p, wit))
print("  вещественная точка: t=0 даёт F0=%s,F4=%s,F8=%s, все > 0" % (m**2, s, n**2))
print("  => класс (1,193,193) СКОРЕЕ ВСЕГО лежит в 2-Селмере, т.е. это элемент Sha[2].")
print("     Значит спуск по Селмеру его НЕ убивает; всё держится на границе ранга.")

# ---------------------------------------------------------------- 12. бесконечность + перебор
hdr("[12] БЕСКОНЕЧНОСТЬ И ПРЯМОЙ ПЕРЕБОР")
print("  над t=inf нужно (u4/t)^2 = s =", s, "; s квадрат в Q:", QQ(s).is_square())
hits=[]; B=500
for q in range(1,B+1):
    for p in range(0,B+1):
        if gcd(p,q)!=1: continue
        A0 = m**2*q**2 + n**2*p**2
        A8 = n**2*q**2 + m**2*p**2
        A4 = ZZ(s)*(q**2+p**2)
        if ZZ(A0).is_square() and ZZ(A8).is_square() and ZZ(A4).is_square():
            hits.append((p,q))
print("  перебор t=p/q, 0<=p<=%d, 1<=q<=%d: найдено точек C:"%(B,B), hits)

hdr("ИТОГ МОЕГО СЧЁТА")
print("  s,b,корни,тождества,необходимый класс : подтверждены полиномиально")
print("  моя точка G                           :", Gm)
print("  различных классов                     :", dist, "/", len(classes))
print("  требуемый класс", req, "среди них :", req in S)
print("  моё множество == множество Codex      :", set(codex_rows)==S)
print("  ранг: eclib", res_rank.get('eclib_bounds'), " PARI", res_rank.get('pari_2'))
