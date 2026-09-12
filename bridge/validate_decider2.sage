# -*- coding: utf-8 -*-
# Вторая, более жёсткая валидация: смещаем выборку к случаям БЕЗ точек
# (домножение на неквадратичные константы), чтобы проверить именно
# отрицательную ветвь решающей процедуры.
load('/home/kep/magicKube/bridge/validate_decider.sage')

print()
print("=== усиленная проверка отрицательной ветви ===")
set_random_seed(777)
agree = 0
disagree = 0
nfalse_conf = 0
ntrue = 0
for trial in range(400):
    p = choice([2, 3, 5, 7])
    ns = 3 if p == 2 else least_quadratic_nonresidue(p)
    c = choice([1, p, ns, p*ns, p*p*ns, -1, -ns])
    while True:
        g = R([ZZ.random_element(-12, 12) for _ in range(7)])
        if g.degree() == 6 and g.is_squarefree():
            break
    f = c*g
    d = decide(f, p)
    if d is None:
        continue
    k = 11 if p == 2 else (8 if p == 3 else 6)
    bt, w = brute_Zp(f, p, k)
    if d == bt:
        agree += 1
        if d is False:
            nfalse_conf += 1
        else:
            ntrue += 1
    else:
        disagree += 1
        print("РАСХОЖДЕНИЕ p=%d c=%s f=%s decide=%s brute=%s %s" % (p, c, f, d, bt, w))
print("УСИЛЕННАЯ ВАЛИДАЦИЯ: согласий %d (из них False подтверждено %d, True %d), расхождений %d"
      % (agree, nfalse_conf, ntrue, disagree))
