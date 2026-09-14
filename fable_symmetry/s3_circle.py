#!/usr/bin/env python3
# Fable, 14.09.2026. Симметричные координаты: четыре точки окружности u^2+v^2=2 / четыре прямоугольных треугольника с гипотенузой 1.
from sage.all import *
import json
out={}
R = PolynomialRing(QQ,'t,t1,t2,t3,t4'); t,t1,t2,t3,t4 = R.gens()
F = R.fraction_field()
# параметризация окружности alpha^2+beta^2=2 через точку (1,1): alpha=(t^2-2t-1)/(t^2+1), beta=(1-2t-t^2)/(t^2+1)
alpha = F((t**2-2*t-1)/(t**2+1)); beta = F((1-2*t-t**2)/(t**2+1))
assert alpha**2+beta**2==2
u = (beta**2-alpha**2)/2
print("u(t) =", u.factor(), "; = 1-alpha^2 = beta^2-1:", u==1-alpha**2 and u==beta**2-1)
assert u==F(4*t*(t**2-1)/(t**2+1)**2)
# катеты треугольника с гипотенузой 1: X=(beta+alpha)/2, Y=(beta-alpha)/2
Xl=(beta+alpha)/2; Yl=(beta-alpha)/2
print("X^2+Y^2 =", Xl**2+Yl**2, "; X =", Xl.factor(), "; Y =", Yl.factor(), "; 4*площадь = 2XY =", (2*Xl*Yl).factor(), "== u:", 2*Xl*Yl==u)
assert Xl**2+Yl**2==1 and 2*Xl*Yl==u
# Пифагоровы (m,n): X = 2mn/(m^2+n^2), Y=(m^2-n^2)/(m^2+n^2) при t = ? Проверим: u = 4mn(m^2-n^2)/(m^2+n^2)^2 при t=n/m
Rm = PolynomialRing(QQ,'m,n'); m,n = Rm.gens(); Fm=Rm.fraction_field()
ut = u.subs({t: Fm(n/m)}) if False else F(4*t*(t**2-1)/(t**2+1)**2)
print("u(n/m) =", factor(4*(n/m)*((n/m)**2-1)/((n/m)**2+1)**2))
# симметрии P^1, сохраняющие u с точностью до знака
for name,g in [('t->-t',-t),('t->1/t',1/t),('t->-1/t',-1/t),('t->(1-t)/(1+t)',(1-t)/(1+t)),('t->(t-1)/(t+1)',(t-1)/(t+1))]:
    ug = F(4*g*(g**2-1)/(g**2+1)**2)
    print(f"  {name}: u -> ", "u" if ug==u else ("-u" if ug==-u else str(ug.factor())))
# вырожденные точки: u=0 (t=0,±1,∞) — клетки совпадают с центром; u=±1 — клетка 0: t^2-2t-1=0 или t^2+2t-1=0 -> t=1±sqrt2, -1±sqrt2: НЕ рациональны
print("u=1 при alpha=0: t^2-2t-1=0 -> t=1±sqrt2 (иррационально) — клетка 0 не достигается рациональной точкой окружности с центром 1; ")
print("  но в проективной записи (клетка 0 при противоположной 2) точка есть: (alpha,beta)=(0,sqrt2) не рациональна! Проверка: 0 и 2 — 2 не квадрат.")
# Действительно: клетки 0 и 2 => 2 должно быть квадратом: невозможно над Q при центре 1. Значит 'клетка 0, противоположная 2' — это центр c=2e^2? нет: центр 1 масштаб. Разберёмся:
print("  Перепись говорит о точках с клеткой 0 и противоположной 2c: при c=1 требует 2=квадрат. Над Q такие точки на V соответствуют центру c с 2c квадратом, т.е. c=2k^2: масштабирование клеток на 2 не сохраняет квадратность. => это точки V, у которых x22=... проверим ниже прямым перебором.")

# Поверхность V в координатах (t1..t4): u3=u1+u2, u4=u1-u2. Магический квадрат Бремнера–Саллоуса (7 квадратов, все 8 линий): не на V; квадрат из переписи z=0: t_i in {0,±1,∞}
U = lambda tt: F(4*tt*(tt**2-1)/(tt**2+1)**2)
eq1 = U(t3)-U(t1)-U(t2); eq2 = U(t4)-U(t1)+U(t2)
print("уравнения V в (P^1)^4: числители степеней", eq1.numerator().degrees(), eq2.numerator().degrees())
# пример над R: b=0.3, c=0.1 -> u = (0.3,0.1,0.4,0.2): точки окружности вещественные: alpha=sqrt(1-u) — вещественная невырожденная точка
us=[QQ(3)/10,QQ(1)/10,QQ(4)/10,QQ(2)/10]
cells=[1]+[1-v for v in us]+[1+v for v in us]
print("вещественная точка: u =", us, " клетки", cells, " все >0 и различны:", all(v>0 for v in cells) and len(set(cells))==9)
out['real_point_cells']=[str(v) for v in cells]

# Квадрат Бремнера–Саллоуса (все 8 линий, 7 квадратов): 373^2 289^2 565^2 / 360721 425^2 23^2 / 205^2 527^2 222121
BS=[[373**2,289**2,565**2],[360721,425**2,23**2],[205**2,527**2,222121]]
Sb=sum(BS[0]); a0=BS[1][1]
assert all(sum(r)==Sb for r in BS) and all(sum(BS[i][j] for i in range(3))==Sb for j in range(3)) and BS[0][0]+BS[1][1]+BS[2][2]==Sb and BS[0][2]+BS[1][1]+BS[2][0]==Sb
b0=a0-BS[0][0]; c0=a0-BS[0][2]
print("Б–С: a =",a0," b =",b0," c =",c0, "; u/a:", [QQ(v)/a0 for v in (b0,c0,b0+c0,b0-c0)])
for name,v in [('b',b0),('c',c0),('b+c',b0+c0),('b-c',b0-c0)]:
    al, be = QQ(a0-v), QQ(a0+v)
    print(f"   пара {name}: {al} {'□' if al.is_square() else '✗'}, {be} {'□' if be.is_square() else '✗'}")
out['BS']=dict(a=a0,b=b0,c=c0)
json.dump(out, open('/home/kep/magicKube/fable_symmetry/s3_circle.json','w'), indent=1)
print("OK s3")
