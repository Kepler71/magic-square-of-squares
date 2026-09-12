# -*- coding: utf-8 -*-
# АТАКА НА РАНГ E для (m,n)=(11,4). Полный 2-спуск (eclib) даёт только 1 <= r <= 3.
# Если r >= 2, аргумент Codex «8 из 8 классов» разваливается.
import sys
from sage.all import *
def flush(): sys.stdout.flush()

m=11; n=4
s=QQ(m^2+n^2)/2; b=s*m^2*n^2
f=[ZZ(-4*b), ZZ(-4*s*m^4), ZZ(-4*s*n^4)]
xv=polygen(QQ,'x'); cub=(xv-f[0])*(xv-f[1])*(xv-f[2])
E=EllipticCurve(QQ,[0,ZZ(cub[2]),0,ZZ(cub[1]),ZZ(cub[0])])
print("E =",E); print("N =",E.conductor(), "=", factor(E.conductor()))
Emin=E.minimal_model(); print("минимальная модель:",Emin)
print("изоморфизм E->Emin:", E.isomorphism_to(Emin))
flush()

print()
print("=== 1. eclib: попытки усилить 2-спуск ===")
from sage.libs.eclib.interface import mwrank_EllipticCurve
for (sl, nx) in [(8,6),(12,8),(16,10),(20,12)]:
    mw=mwrank_EllipticCurve(list(E.a_invariants())); mw.set_verbose(0)
    try:
        mw.two_descent(second_limit=sl, n_aux=nx)
    except Exception as ex:
        print("  two_descent(second_limit=%d) ошибка: %s"%(sl,ex))
    print("  second_limit=%2d n_aux=%2d -> rank=%s bound=%s certain=%s selmer=%s"
          %(sl,nx,mw.rank(),mw.rank_bound(),mw.certain(),mw.selmer_rank()))
    flush()
    if mw.certain(): break

print()
print("=== 2. PARI ellrank (независимый 2-спуск) ===")
try:
    er = pari(E).ellrank()
    print("  ellrank =", er)
    print("  интерпретация: [нижняя, верхняя, s, точки]")
except Exception as ex:
    print("  ошибка:", ex)
flush()
try:
    er2 = pari(Emin).ellrank()
    print("  ellrank(мин.модель) =", er2)
except Exception as ex:
    print("  ошибка:", ex)
flush()

print()
print("=== 3. Аналитический ранг (через него — Колывагин/Гросс–Загир) ===")
for alg in ['pari','rubinstein','sympow']:
    try:
        ar = E.analytic_rank(algorithm=alg)
        print("  analytic_rank(%s) = %s" % (alg, ar))
    except Exception as ex:
        print("  analytic_rank(%s): ошибка %s" % (alg, ex))
    flush()
try:
    ub = E.analytic_rank_upper_bound()
    print("  analytic_rank_upper_bound() =", ub, " (строгая верхняя граница)")
except Exception as ex:
    print("  analytic_rank_upper_bound: ошибка", ex)
flush()
try:
    print("  знак функционального уравнения w =", E.root_number())
except Exception as ex:
    print("  root_number: ошибка", ex)

print()
print("=== 4. Поиск дополнительных независимых точек (ключ: r>=2 ломает вывод) ===")
found=[]
for hb in [12,16,20,24]:
    try:
        pts = E.point_search(hb, verbose=False)
        print("  point_search(%d): %d точек" % (hb, len(pts)))
        for P in pts:
            if P.order()==Infinity: found.append(P)
        flush()
    except Exception as ex:
        print("  point_search(%d): ошибка %s" % (hb,ex)); flush()
    if found:
        try:
            sat,idx,reg = E.saturation(found)
            print("    насыщение: %d образующих, индекс %s, регулятор %s" % (len(sat),idx,reg))
            hm = E.height_pairing_matrix(sat)
            print("    ранг решётки найденных точек =", hm.rank())
            for P in sat: print("      ", P, " h=", RR(P.height()))
        except Exception as ex:
            print("    saturation: ошибка", ex)
        flush()

print()
print("=== 5. Sha и BSD-эвристика ===")
try:
    print("  Sha.an_numerical() =", E.sha().an_numerical())
except Exception as ex:
    print("  Sha: ошибка", ex)
try:
    print("  Sha.an() =", E.sha().an())
except Exception as ex:
    print("  Sha.an(): ошибка", ex)
flush()

print()
print("=== 6. Лежит ли ТРЕБУЕМЫЙ класс (1,274,274) в 2-Селмере? ===")
print("    Если да и если ранг окажется 3, то образ = Селмер и класс ДОСТИГАЕТСЯ.")
print("    Локальный тест: существует ли x в Q_p c (x-f1)/1, (x-f2)/274, (x-f3)/274 квадратами.")
d=[ZZ(1),ZZ(274),ZZ(274)]
def is_sq_Qp(a,p,prec_ok=True):
    a=QQ(a)
    if a==0: return True
    v=a.valuation(p)
    if v%2: return False
    u=a/p**v
    if p==2:
        return (ZZ(u.numerator()*u.denominator()) % 8) == 1
    num=ZZ(u.numerator()); den=ZZ(u.denominator())
    return kronecker(num*den % p, p)==1 if (num*den)%p!=0 else False

def local_solvable(p, kmax=3, nmax=None):
    """поиск x в Q_p (через рациональные приближения) с нужными квадратичными классами"""
    if nmax is None: nmax = 6 if p>2 else 9
    # x = p^a * u, перебираем a и u по модулю p^nmax
    for a in range(-4, 9):
        pa = QQ(p)**a
        rng = p**nmax
        for u in range(1, rng):
            if u % p == 0: continue
            for sg in ([1,-1] if True else [1]):
                x = sg*pa*u
                ok=True
                for i in range(3):
                    val = (x - f[i])/d[i]
                    if val == 0:
                        ok=False; break
                    # достаточная точность: сравниваем с точным рациональным x
                    if not is_sq_Qp(val,p):
                        ok=False; break
                if ok:
                    return (True, x)
    return (False, None)

bad = sorted(set(E.conductor().prime_factors()+[2]))
print("  плохие простые:", bad)
for p in bad:
    res, wit = local_solvable(p)
    print("    p=%-4s локально достижим ? %s   свидетель x=%s" % (p,res,wit))
    flush()
# вещественное место
print("    R: нужно x >= max(f) =", max(f), "и все три положительны -> да (неограниченная компонента)")
print("    ВНИМАНИЕ: этот локальный поиск — перебор, отрицательный ответ не окончателен.")
