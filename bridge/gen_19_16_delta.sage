# -*- coding: utf-8 -*-
# G1 (19,16): точные классы delta, торсион, и локальная разрешимость требуемого класса.
# Всё — точной рациональной арифметикой.
import sys, itertools

m = 19; n = 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
e = [-b, -s*m^4, -s*n^4]
print("s=",s," b=",b," e=",e)

R.<X> = QQ[]
cub = prod(X - ei for ei in e)
E0 = EllipticCurve([0, cub[2], 0, cub[1], cub[0]])
Emin = EllipticCurve([0,0,0,-1613212693425900,7222574158576522830000])
print("E0=",E0); print("Emin=",Emin)

# точная проверка изоморфизма E0 -> Emin: x_min = 4X + 118576294, y_min = 8V
u = QQ(1)/2          # X = u^2*x + r
r = QQ(-59288147)/2
# проверяем тождество полиномов
xx = polygen(QQ,'xx')
lhs = (xx^3 - 1613212693425900*xx + 7222574158576522830000)   # = y_min^2
# V^2 = cub(X), X = u^2 xx + r, V = u^3 y_min => u^6 y^2 = cub(u^2 xx + r)
chk = cub(u^2*xx + r) - u^6*lhs
print("проверка изоморфизма (должно быть 0):", chk)
assert chk == 0
print("u^2 =", u^2, " квадрат в Q*: ", (u^2).is_square())

# корни Emin
Emin_roots = sorted((xx^3 - 1613212693425900*xx + 7222574158576522830000).roots(multiplicities=False))
print("корни Emin (sorted):", Emin_roots)
# согласованный порядок: E_i = (e_i - r)/u^2 = 4*e_i + 118576294
Ecoord = [(ei - r)/u^2 for ei in e]
print("корни Emin в порядке e1,e2,e3:", Ecoord)
assert sorted(Ecoord) == Emin_roots

def sqclass(q):
    q = QQ(q)
    assert q != 0
    return q.squarefree_part()

# ТРЕБУЕМЫЙ КЛАСС
req = (sqclass(1), sqclass(s), sqclass(s))
print("ТРЕБУЕМЫЙ КЛАСС delta = (1, s, s) =", req)

# --- delta для точки на Emin (использует Ecoord = 4e_i+118576294, масштаб 4 = квадрат) ---
E1c, E2c, E3c = Ecoord
def delta_point(P):
    if P.is_zero():
        return (1,1,1)
    x = P[0]; y = P[1]
    d = []
    Ec = [E1c,E2c,E3c]
    for i in range(3):
        v = x - Ec[i]
        if v == 0:
            j,k = [t for t in range(3) if t != i]
            v = (Ec[i]-Ec[j])*(Ec[i]-Ec[k])
        d.append(sqclass(v))
    return tuple(d)

tors = Emin.torsion_points()
print("\nторсион Emin:", [ (P[0],P[1]) if not P.is_zero() else 'O' for P in tors])
tors_cls = {}
for P in tors:
    c = delta_point(P)
    tors_cls[c] = P
    print("  delta(",("O" if P.is_zero() else (P[0],P[1])),") =", c)
print("различных торсионных классов:", len(tors_cls))
print("требуемый класс среди торсионных?", req in tors_cls)
sys.stdout.flush()

# ---------------------------------------------------------------
# ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ требуемого класса (2-Селмер, полная 2-спуск-модель)
# Однородное пространство для (d1,d2,d3):
#   d1 z1^2 - d2 z2^2 = E2 - E1
#   d1 z1^2 - d3 z3^2 = E3 - E1
# (стандартная модель полного 2-спуска, корни Ecoord)
# ---------------------------------------------------------------
d1,d2,d3 = [QQ(t) for t in req]
A = E2c - E1c
B = E3c - E1c
print("\nмодель однородного пространства для требуемого класса:")
print("  d=(%s,%s,%s)"%(d1,d2,d3))
print("  %s z1^2 - %s z2^2 = %s" % (d1,d2,A))
print("  %s z1^2 - %s z3^2 = %s" % (d1,d3,B))

def local_sol(d1,d2,d3,A,B,p,prec=None):
    """грубая проверка p-адической разрешимости перебором по Z_p/p^k"""
    if p == 0:
        # вещественно: нужны z1,z2,z3 вещественные
        # d1 z1^2 - d2 z2^2 = A ;  d1 z1^2 - d3 z3^2 = B
        # берём z1^2 = w >= 0 ; нужно (d1 w - A)/d2 >= 0 и (d1 w - B)/d3 >= 0
        import numpy
        for w in [QQ(k)/16 for k in range(0,100000)]:
            if (d1*w-A)/d2 >= 0 and (d1*w-B)/d3 >= 0:
                return True
        return False
    k = prec or (6 if p>2 else 12)
    Rp = Zp(p, prec=k+8, type='capped-rel')
    # ищем решение в Q_p: перебираем z1 по представителям mod p^k с валюациями
    # используем Hensel через проверку квадратов
    for v1 in range(-4,5):
        for u1 in range(0, p^k):
            if u1 % p == 0 and v1 > -4:
                continue
            z1sq = QQ(u1)*QQ(p)^(2*v1)
            if z1sq == 0: continue
            t2 = (d1*z1sq - A)/d2
            t3 = (d1*z1sq - B)/d3
            if t2 == 0 or t3 == 0:
                if t2 == 0 and t3 == 0: return True
                other = t3 if t2==0 else t2
                if Qp(p,k+8)(other).is_square(): return True
                continue
            if Qp(p,k+8)(t2).is_square() and Qp(p,k+8)(t3).is_square():
                return True
    return None   # не нашли — неубедительно

print("\n(грубая локальная проверка выполняется отдельно ниже)")
