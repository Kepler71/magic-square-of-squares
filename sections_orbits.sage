# Галуа-орбиты целых сечений: минимальные ассоциированные простые идеала системы.
import time
Rs.<s> = QQ[]
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
B.<a,b,c,d,e,f,g> = QQ[]
Bs.<t> = B[]
X = a*t^2 + b*t + c; Y = d*t^3 + e*t^2 + f*t + g
I = B.ideal((Y^2 - (X^3 + Bs(A4)(t)*X + Bs(A6)(t))).coefficients())
t0 = time.time()
J = I.radical()
print("distinct integral sections over Qbar:", J.vector_space_dimension(), f"({time.time()-t0:.1f}s)")
mp = J.minimal_associated_primes()
print("Galois orbits:", len(mp), " sizes:", sorted(Pj.vector_space_dimension() for Pj in mp))
save(mp, 'sections_orbits.sobj')
for Pj in sorted(mp, key=lambda Pj: Pj.vector_space_dimension()):
    dj = Pj.vector_space_dimension()
    # поле определения: минимальный многочлен порождающей линейной формы
    for lin in [a, c, b + 3*c, d + 5*e + 7*g, a + 2*b + 3*c + 5*d + 7*e + 11*f + 13*g]:
        Z.<w> = QQ[]
        Jl = (Pj + B.ideal(B('w0') if False else [])).gens()
        # исключаем всё, кроме lin: вводим новую переменную
        C.<a2,b2,c2,d2,e2,f2,g2,w2> = QQ[]
        phi = B.hom([a2,b2,c2,d2,e2,f2,g2], C)
        K = C.ideal([phi(h) for h in Pj.gens()] + [w2 - phi(lin)])
        m = K.elimination_ideal([a2,b2,c2,d2,e2,f2,g2]).gens()[0]
        if m.degree() == dj:
            break
    mpol = Z(m.univariate_polynomial()(w)).monic()
    Kf = NumberField(mpol, 'z') if dj > 1 else QQ
    info = f"disc {Kf.discriminant().factor()}, Galois closure deg {Kf.galois_closure('z2').degree() if dj > 1 else 1}" if dj > 1 else "rational"
    sol_a = [h for h in Pj.gens()]
    print(f"--- orbit size {dj}: field {info}")
    if dj <= 2:
        for sol in Pj.variety(QQbar if dj > 1 else QQ):
            print("     X =", sol[a], "s^2 +", sol[b], "s +", sol[c], "   Y-coeffs:", sol[d], sol[e], sol[f], sol[g])
    else:
        print("     min poly of generator:", mpol)
