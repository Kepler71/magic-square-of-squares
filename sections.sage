# Задача C: все "целые" сечения (не пересекающие O) лучевой поверхности
#   Y^2 = X^3 - 432 P(s) X + 3456 (s-6)(3s-2) Q(s)
# Для рациональной эллиптической поверхности такие сечения: deg X <= 2, deg Y <= 3.
# Подстановка -> 7 уравнений на 7 неизвестных; нульмерная система, решаем Грёбнером над Q.
import time
Rs.<s> = QQ[]
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4 = -432*P
A6 = 3456*(s-6)*(3*s-2)*Q
print("deg a4, a6:", A4.degree(), A6.degree())
cub = lambda X: X^3 + A4*X + A6

# контроль известных сечений из записки
T  = (36*s^2 - 240*s + 144, 0)
S1 = (72*s^2 - 96*s + 288, 432*(s+2)*(s^2-12))
print("T on surface:", T[1]^2 == cub(T[0]))
print("S1 on surface:", S1[1]^2 == cub(S1[0]))
Fs = Rs.fraction_field()
TS = (Fs(12*(3*s^4-8*s^3-248*s^2-32*s+1584)/(s+2)^2), Fs(27648*(s^2-12)*(s^2-8)/(s+2)^3))
print("T+S1 on surface:", TS[1]^2 == cub(TS[0]))
E = EllipticCurve(Fs, [A4, A6])
print("T+S1 == T + S1 in E(Q(s)):", E(T) + E(S1) == E(TS))
print("disc factor:", E.discriminant().factor())

# система на целые сечения
B.<a,b,c,d,e,f,g> = QQ[]
Bs.<t> = B[]
X = a*t^2 + b*t + c
Y = d*t^3 + e*t^2 + f*t + g
eq = Y^2 - (X^3 + Bs(A4)(t)*X + Bs(A6)(t))
I = B.ideal(eq.coefficients())
t0 = time.time()
print("\ndim:", I.dimension(), " #solutions over Qbar (with mult.):", I.vector_space_dimension(), f"({time.time()-t0:.1f}s)")
# поле определения: исключаем всё, кроме a и c (координата X определяет сечение с точностью до знака Y)
Ja = I.elimination_ideal([b, c, d, e, f, g])
pa = Ja.gens()[0].univariate_polynomial()
print("elim poly in a: deg", pa.degree())
for fac, m in pa.factor():
    print("   ", fac, "  mult", m)
print("\nrational points of the system:")
for sol in I.variety(QQ):
    print("   X =", sol[a]*s^2+sol[b]*s+sol[c], "  Y =", sol[d]*s^3+sol[e]*s^2+sol[f]*s+sol[g])
