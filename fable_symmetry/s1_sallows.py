#!/usr/bin/env python3
# Fable, 14.09.2026. Квадрат Саллоуса (7 линий, 9 квадратов): проверки и структура.
from sage.all import *
import itertools, json

M = [[127,46,58],[2,113,94],[74,82,97]]
C = [[x*x for x in r] for r in M]
S = sum(C[0])
out = {}
print("клетки:", C)
rows = [sum(r) for r in C]; cols=[sum(C[i][j] for i in range(3)) for j in range(3)]
d1 = C[0][0]+C[1][1]+C[2][2]; d2 = C[0][2]+C[1][1]+C[2][0]
print("строки", rows, "столбцы", cols, "диаг", d1, "антидиаг", d2)
assert rows==[S]*3 and cols==[S]*3 and d2==S and d1!=S
print("S =", S, "= 147^2:", S==147**2, "; центр 113^2, S/3 =", QQ(S)/3, "квадрат?", QQ(S/3).is_square())
print("центральная симметрия: 127^2+97^2 =", 127**2+97**2, " 2*113^2 =", 2*113**2)
assert len(set(sum(C,[])))==9
out['sums_ok']=True

# 1. Знаковая ортогональность: существует ли выбор знаков eps_ij, при котором A=(eps_ij*M_ij) конформна (A A^T = S I)?
sols=[]
for eps in itertools.product([1,-1],repeat=9):
    A = matrix(ZZ,3,3,[eps[3*i+j]*M[i][j] for i in range(3) for j in range(3)])
    if A*A.transpose()==S*identity_matrix(3):
        sols.append(A)
print("знаковых конформных вариантов:", len(sols))
A = sols[0]; print(A); print("A A^T =", (A*A.transpose()).diagonal(), "det =", A.det(), "S^(3/2) =", 147**3)
out['signed_orthogonal_variants']=len(sols)
# все варианты = один с точностью до знаков строк и столбцов (2^3*2^3/2=32)?
assert len(sols)==32

# 2. Разложение 7-линейного квадрата: C = Lucas(a,b,c) + lam*I
lam = QQ(3*C[1][1]-S)/2
a = C[1][1]-lam
Mag = [[C[i][j]-(lam if i==j else 0) for j in range(3)] for i in range(3)]
print("lambda =", lam, " a =", a, " магическая часть:", Mag)
Sm = sum(Mag[0])
assert all(sum(r)==Sm for r in Mag) and all(sum(Mag[i][j] for i in range(3))==Sm for j in range(3))
assert Mag[0][0]+Mag[1][1]+Mag[2][2]==Sm and Mag[0][2]+Mag[1][1]+Mag[2][0]==Sm and Sm==3*a
b = a-Mag[0][0]; c = a-Mag[0][2]
print("b =", b, " c =", c, " проверка Лукаса:", Mag==[[a-b,a+b+c,a-c],[a+b-c,a,a-b+c],[a+c,a-b-c,a+b]])
assert Mag==[[a-b,a+b+c,a-c],[a+b-c,a,a-b+c],[a+c,a-b-c,a+b]]
out.update(dict(lam=int(lam),a=int(a),b=int(b),c=int(c)))

# 3. Четыре пары через центр: три с серединой a, одна (главная диагональ) с серединой a+lam
pairs = {'диаг (b)':(C[0][0],C[2][2]), 'антидиаг (c)':(C[0][2],C[2][0]), 'столбец (b+c)':(C[0][1],C[2][1]), 'строка (b-c)':(C[1][0],C[1][2])}
for k,(x,y) in pairs.items(): print("  пара", k, ": середина", QQ(x+y)/2, " полуразность", QQ(y-x)/2)

# 4. b = -3360 = -16*210; 210 конгруэнтно (треугольник 20,21,29); кривая E_b: y^2=x^3-b^2 x имеет ранг>=1
n = 210
E = EllipticCurve([-n**2,0])
print("E_210: y^2=x^3-210^2 x, ранг (PARI, точно?)", E.rank(), " кручение", E.torsion_order())
# точка на E_b из тройки квадратов a+lam-b=127^2, a+lam=113^2, a+lam+b=97^2 (AP с разностью b):
# стандарт: для AP alpha^2, gamma^2, beta^2 с разностью n: x = gamma^2 (после масштаба) -> точка (x, y) c y^2 = x^3 - n^2 x, где x=gamma^2, n = разность
g2 = 113**2; nn = 3360
E2 = EllipticCurve([-nn**2,0])
P = E2.lift_x(g2)
print("точка на E_3360 с x = 113^2:", P, " порядок:", P.order())
# новые 7-линейные квадраты с теми же шестью клетками: x(2P), x(3P) должны давать AP квадратов gamma^2-3360, gamma^2, gamma^2+3360
def ap_from_point(Q):
    x = Q[0]
    assert x.is_square() and (x-nn).is_square() and (x+nn).is_square()
    return x
new=[]
for k in [2,3]:
    Q = k*P
    x = ap_from_point(Q)
    # квадрат: клетки a±c, a±(b+c), a±(b-c) масштабируем? Нет: нужны та же шестёрка -> центр a+lam' = x требует того же масштаба.
    # x рационально; клетки: (1,1)=x-b, (2,2)=x, (3,3)=x+b, остальные как у Саллоуса.
    sq = [[x-b, C[0][1], C[0][2]],[C[1][0], x, C[1][2]],[C[2][0], C[2][1], x+b]]
    ok = all(QQ(v).is_square() for r in sq for v in r) and all(sum(r)==sum(sq[0]) for r in sq) and all(sum(sq[i][j] for i in range(3))==sum(sq[0]) for j in range(3)) and sq[0][2]+sq[1][1]+sq[2][0]==sum(sq[0])
    print(f"  {k}P: x = {x} = ({sqrt(QQ(x))})^2 ; новый 7-линейный квадрат из 9 квадратов, различных: {ok and len(set(QQ(v) for r in sq for v in r))==9}")
    print("     корни:", [[sqrt(QQ(v)) for v in r] for r in sq])
    new.append([[str(QQ(v)) for v in r] for r in sq])
    assert ok
out['new_seven_line_squares']=new
json.dump(out, open('/home/kep/magicKube/fable_symmetry/s1_sallows.json','w'), ensure_ascii=False, indent=1)
print("OK s1")
