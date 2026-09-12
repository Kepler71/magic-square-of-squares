# Q-кривая: E1 изогенна своей галуа-сопряжённой над k-баром. Достаточное условие: j(E1) в Q
# (тогда E1 — твист кривой над Q, сопряжённая — тоже твист, и они изоморфны над расширением).
# Тогда J2 = Res_{k/Q} E1 имеет вещественное умножение и rk NS(J2) = 2 (условие квадратичного Шаботи: rank < g + rho - 1 = 3).
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
data = [('(17,7,13) s=1/5', 109, 229, 65), ('(7,1,5)', 13, 37, 5), ('(23,7,17)', 169, 409, 34),
        ('(31,17,25)', 457, 793, 7), ('(41,1,29)', 421, 1261, 29), ('(47,23,37)', 949, 1789, 37),
        ('(49,31,41)', 1321, 2041, 41), ('(73,17,53)', 1549, 4069, 53), ('(71,49,61)', 3061, 4381, 671),
        ('(89,23,65)', 2377, 6073, 910), ('(79,47,65)', 3217, 5233, 455)]
for name, A, C, D in data:
    try:
        k, rts = E1_roots(A, C, D)
        E = EllipticCurve(k, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)])
        j = E.j_invariant()
        inQ = j in QQ
        sig = k.embeddings(k)[1] if len(k.embeddings(k)) > 1 else None
        Es = EllipticCurve(k, [sig(c) for c in E.ainvs()]) if sig else None
        iso = Es.is_isomorphic(E) if Es else None
        print(f"{name}: j в Q: {inQ}; сопряжённая изоморфна над k: {iso}")
    except Exception as ex:
        print(f"{name}: ошибка {type(ex).__name__}: {str(ex)[:60]}")
