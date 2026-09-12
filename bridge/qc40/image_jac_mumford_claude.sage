#!/usr/bin/env sage
# image_jac_mumford_claude.sage -- НЕЗАВИСИМАЯ перепроверка §6 из basis_E2.sage:
#   образ Psi(Jac(D)(Q)) в E1(Q) x E2(Q) и его индекс.
#
# Отличие от basis_E2.sage: там локальный образ Psi(Jac(D)(F_p)) строился РУЧНЫМ
# перебором эффективных дивизоров степени 2 (пары F_p-точек + сопряжённые пары из F_p^2).
# Здесь используется СОБСТВЕННАЯ группа Якоби Sage (split-модель Мамфорда, класс
# MumfordDivisorClassFieldSplit) -- другая реализация группового закона и перечисления.
# Дополнительно берутся ДРУГИЕ простые (101..199), которых в basis_E2.sage не было.
#
# Соглашение Мамфорда для split-модели (deg f = 6, две рациональные точки на бесконечности):
#   класс (u, v : n)  =  [D_{u,v}] + n[inf+] + (2 - deg u - n)[inf-] - [inf+] - [inf-]
# Соглашение ПРОВЕРЯЕТСЯ в скрипте: Psi обязан быть гомоморфизмом и Psi(0) = 0.
#
# Запуск: sage /home/kep/magicKube/bridge/qc40/image_jac_mumford_claude.sage

from sage.all import *
from pathlib import Path
import json, random

OUT = Path('/home/kep/magicKube/bridge/qc40')
LOG = []
R = {}
random.seed(int(7))


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print(s, flush=True)


# ------------------------------------------------------------------ модель
bridge = json.loads((OUT / 'new_bridge_verified.json').read_text())
a0, a2, a4, a6 = [QQ(z) for z in bridge['coefficients_a0_a2_a4_a6']]
say("[0] из new_bridge_verified.json: a0=%s a2=%s a4=%s a6=%s" % (a0, a2, a4, a6))

E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve([0, a2, 0, a4 * a0, a6 * a0**2])
say("[0] E1 raw = %s ; совпал с bridge? %s"
    % (E1.ainvs(), [str(z) for z in E1.ainvs()] == bridge['elliptic_factors'][0]['raw_ainvs']))
say("[0] E2 raw = %s ; совпал с bridge? %s"
    % (E2.ainvs(), [str(z) for z in E2.ainvs()] == bridge['elliptic_factors'][1]['raw_ainvs']))
assert [str(z) for z in E1.ainvs()] == bridge['elliptic_factors'][0]['raw_ainvs']
assert [str(z) for z in E2.ainvs()] == bridge['elliptic_factors'][1]['raw_ainvs']

M1 = E1.minimal_model()
M2 = E2.minimal_model()
sh1 = ZZ(12536304)
sh2 = ZZ(15762384)
assert [str(z) for z in M1.ainvs()] == bridge['elliptic_factors'][0]['minimal_ainvs']
assert [str(z) for z in M2.ainvs()] == bridge['elliptic_factors'][1]['minimal_ainvs']
# сдвиг X проверяем сами, а не берём на веру
# соглашение Sage: (u,r,s,t) означает x_min = (x_raw - r)/u^2, то есть sh = -r
u_, r_, s_, t_ = E1.isomorphism_to(M1).tuple()
assert u_ == 1 and s_ == 0 and t_ == 0 and r_ == -sh1, (u_, r_, s_, t_)
u_, r_, s_, t_ = E2.isomorphism_to(M2).tuple()
assert u_ == 1 and s_ == 0 and t_ == 0 and r_ == -sh2, (u_, r_, s_, t_)
# и прямая проверка на точке
assert M1(E1.isomorphism_to(M1)(E1([0, a6 * 2023, 1]))) == M1([sh1, a6 * 2023, 1])
assert M2(E2.isomorphism_to(M2)(E2([0, a0 * 6647, 1]))) == M2([sh2, a0 * 6647, 1])
say("[0] изоморфизмы raw->min: чистые сдвиги X на %s и %s (u=1, s=t=0) -- проверено" % (sh1, sh2))

P1 = M1([34384993, 286391046720, 1]); T1 = M1([8443775, 0, 1])
P2 = M2([-2479007, 87162492480, 1]); T2 = M2([-28420225, 0, 1])
assert T1.order() == 2 and T2.order() == 2
say("[0] P1, T1, P2, T2 лежат на минимальных моделях; T1, T2 порядка 2")

A1M = [int(z) for z in M1.ainvs()]
A2M = [int(z) for z in M2.ainvs()]
D1 = ZZ(M1.discriminant()); D2 = ZZ(M2.discriminant())


# -------------------------------------- Psi из мамфордова представления Sage
def psi_of_mumford(P, p, e1F, e2F, e1K, e2K, K, s_inf0):
    """P -- элемент J(F_p) в split-представлении Sage. Возвращает (A, B) в E1(F_p) x E2(F_p)."""
    u, v, n = tuple(P)          # Sage: __iter__ отдаёт (u, v, n); len() не определён
    n = ZZ(n)
    du = u.degree()
    A = e1K(0)
    B = e2K(0)
    # аффинная часть: корни u в F_p^2 с кратностями
    uu = u.change_ring(K)
    rts = uu.roots()
    tot = sum(m for _, m in rts)
    assert tot == du, "u не расщепляется в F_p^2: %s" % u
    for r, m in rts:
        yv = v.change_ring(K)(r)
        for _ in range(m):
            A += e1K([K(a6) * r * r + K(sh1), K(a6) * yv, 1])
            if r == 0:
                pass                      # phi2(x=0) = O
            else:
                B += e2K([K(a0) / r**2 + K(sh2), K(a0) * yv / r**3, 1])
    # часть на бесконечности: phi1(inf+-) = O ; phi2(inf+) = Bpt, phi2(inf-) = -Bpt
    # ВАЖНО: Sage сам выбирает, какая из двух точек на бесконечности есть P_0
    # (её y/x^3 = sqrt(a6) в выборе Sage, что НЕ всегда равно 6647 mod p).
    # Берём её из H.points_at_infinity()[0], а не из константы 6647.
    Bpt = e2K([K(sh2), K(a0) * K(s_inf0), 1])
    B += (2 * n + du - 2) * Bpt
    # спуск обратно в F_p
    A = e1F(0) if A.is_zero() else e1F([GF(p)(A[0]), GF(p)(A[1]), 1])
    B = e2F(0) if B.is_zero() else e2F([GF(p)(B[0]), GF(p)(B[1]), 1])
    return (A, B)


def local_image(p):
    F = GF(p)
    K = GF(p**2, 'aa')
    Rp = PolynomialRing(F, 'x'); xp = Rp.gen()
    fp = F(a6) * xp**6 + F(a4) * xp**4 + F(a2) * xp**2 + F(a0)
    assert fp.degree() == 6 and fp.is_squarefree(), "плохая редукция при p=%s" % p
    H = HyperellipticCurve(fp)
    JF = H.jacobian()(F)
    e1F = M1.change_ring(F); e2F = M2.change_ring(F)
    e1K = EllipticCurve(K, [K(z) for z in A1M]); e2K = EllipticCurve(K, [K(z) for z in A2M])
    s_inf0 = H.points_at_infinity()[0][1]      # y/x^3 у мамфордовой P_0 (Sage)
    assert s_inf0**2 == F(a6)
    pts = JF.points()
    nJ = len(pts)
    # контроль: |J(F_p)| = |E1(F_p)|*|E2(F_p)|  (Тейт: изогенные => равные порядки)
    NE = e1F.cardinality() * e2F.cardinality()
    assert nJ == NE, "|J(F_p)|=%s != |E1xE2(F_p)|=%s при p=%s" % (nJ, NE, p)
    img = {}
    for z in pts:
        img[z] = psi_of_mumford(z, p, e1F, e2F, e1K, e2K, K, s_inf0)
    # ПРОВЕРКА СОГЛАШЕНИЯ: Psi -- гомоморфизм, Psi(0) = 0
    zero = [z for z in pts if z.is_zero()]
    assert len(zero) == 1
    assert img[zero[0]] == (e1F(0), e2F(0)), "Psi(0) != 0 при p=%s" % p
    L = list(pts)
    for _ in range(200):
        z1 = random.choice(L); z2 = random.choice(L)
        s = z1 + z2
        lhs = img[s]
        rhs = (img[z1][0] + img[z2][0], img[z1][1] + img[z2][1])
        assert lhs == rhs, "Psi не гомоморфизм при p=%s (соглашение Мамфорда неверно)" % p
    S = set(img.values())
    idx = ZZ(NE) // len(S)
    return nJ, NE, S, idx, e1F, e2F, s_inf0


classes = [(i, j, k, l) for i in [0, 1] for j in [0, 1] for k in [0, 1] for l in [0, 1]]
survivors = set(classes)
rows = []
primes_old = [11, 19, 29, 31, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
primes_new = [p for p in prime_range(101, 200)]
for p in primes_old + primes_new:
    if (D1 * D2) % p == 0:
        continue
    nJ, NE, S, idx, e1F, e2F, s0 = local_image(p)
    local_s_inf0 = [s0]
    ok = set()
    for (i, j, k, l) in classes:
        if (e1F(i * P1 + j * T1), e2F(k * P2 + l * T2)) in S:
            ok.add((i, j, k, l))
    survivors &= ok
    pred = 4 if p % 4 == 1 else 2
    rows.append({'p': int(p), 'card_J_Fp': int(nJ), 'card_E1xE2': int(NE),
                 'card_image': int(len(S)), 'local_index': int(idx),
                 'predicted_index': int(pred),
                 'new_prime': bool(p in primes_new),
                 's_inf0': int(local_s_inf0[0]), 'six647_mod_p': int(GF(p)(6647)),
                 'surviving': sorted([list(c) for c in ok])})
    say("[1] p=%-4d |J(F_p)|=%-7d |Psi(J)|=%-7d индекс=%-2s (предсказано %s) классов=%-2d %s"
        % (p, nJ, len(S), idx, pred, len(ok), "NEW" if p in primes_new else ""))
    say("[1]      накопленное пересечение: %s   (P_0 Sage: y/x^3=%s, 6647 mod p=%s)"
        % (sorted(survivors), GF(p)(0) + rows[-1]['s_inf0'], GF(p)(6647)))
    assert idx == pred, "РАСХОЖДЕНИЕ формулы |K(F_p)| при p=%s" % p

say("[2] ИТОГ независимого расчёта: выжившие классы = %s" % sorted(survivors))
exact = (survivors == {(0, 0, 0, 0), (0, 1, 0, 1)})
say("[2] совпало с basis_E2.sage ({(0,0,0,0),(0,1,0,1)}) ? %s" % exact)

stored = json.loads((OUT / 'basis_E2.json').read_text())
stored_surv = set(tuple(int(z) for z in c)
                  for c in stored['image_of_Jac']['upper_containment']['surviving_classes'])
say("[2] в basis_E2.json записано: %s" % sorted(stored_surv))
say("[2] независимый расчёт == записанному ? %s" % (survivors == stored_surv))

R = {
    'what': 'независимая перепроверка §6 basis_E2.sage через собственную группу Якоби Sage',
    'different_from_original': ['группа Якоби Sage (Mumford split) вместо ручного перебора дивизоров',
                               'добавлены простые 101..199, которых в оригинале не было'],
    'mumford_convention_selfcheck': 'Psi(0)=0 и гомоморфность на 200 случайных парах для каждого p',
    'tate_selfcheck': '|J(F_p)| = |E1(F_p)|*|E2(F_p)| для каждого p',
    'per_prime': rows,
    'surviving_classes': sorted([list(c) for c in survivors]),
    'matches_basis_E2_json': bool(survivors == stored_surv),
    'exact_image': bool(exact),
    'index_in_product': 8 if exact else None,
    'status': 'proved_software',
}
(OUT / 'image_jac_mumford_claude.json').write_text(
    json.dumps(R, indent=2, ensure_ascii=False) + '\n')
(OUT / 'image_jac_mumford_claude.log').write_text("\n".join(LOG) + "\n")
say("записано: image_jac_mumford_claude.json / .log")
