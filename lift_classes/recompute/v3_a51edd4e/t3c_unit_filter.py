# -*- coding: utf-8 -*-
"""t3c: СТРОГИЙ фильтр по единичным частям при простых p = 3 mod 4, p >= 7 (своё уточнение теоремы Codex).

Лемма (доказательство -- в отчёте; коротко):
  (a) p = 3 mod 4: у z нет полюса. Пусть lam из {r,s} -- p-единица (r,s взаимно просты).
      При v_p(z) = -n < 0: [1 +- lam z] = [+- lam z] = (n, chi(+-lam u)); по матрице S/U классы 1+lam z и 1-lam z
      равны => chi(-1)=1, противоречие.
  (b) z в Z_p: для каждой пары противоположных клеток 1 +- mu z (сумма = 2) хотя бы одна -- единица,
      её класс = (0, chi(вычет)). Поэтому классы U (mu=r,s), S (mu=r+s), S+U (mu=r-s) определяются z0 = z mod p,
      и восемь произведений возможны лишь при z0, для которых эти значения согласованы:
      chi-классы двух членов пары (если оба единицы) равны, c_r = c_s, c_{r+s} c_{r-s} = c_r.
  => [T]_p = c_{r+s}(z0), [L]_p = c_{r-s}(z0) для некоторого согласованного z0 (z0 = 0 даёт (1,1)).
  => для глобальных d_T, d_L (они p-единицы): (chi_p(d_T), chi_p(d_L)) лежит в I_p = {(c_{r+s}(z0), c_{r-s}(z0))}.
Это ВЕРХНЯЯ оценка образа; для исключения кандидата её достаточно.

Скрипт: (1) для трёх наклонов -- подробности; (2) для всех наклонов s <= SMAX с нетривиальными списками Codex --
  уточнённые списки по простым 7 <= p <= PMAX, p = 3 mod 4; (3) сверка: выборочный точный образ (t3b-метод,
  реальные рациональные z) содержится в I_p, а также I_p реализуется (нижняя оценка) -- на контрольных наклонах.
"""
import json, sys, time, random
from math import gcd
from t1_lists import admissible_A

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 300

def primes_upto(n):
    sv = bytearray([1]) * (n + 1); sv[0] = sv[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(n + 1) if sv[i]]

P3 = [p for p in primes_upto(PMAX) if p % 4 == 3 and p >= 7]

def chi(x, p):
    x %= p
    if x == 0: return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1

def image_upper(r, s, p):
    """I_p как множество пар (+-1,+-1) = (chi-класс [T], chi-класс [L])"""
    img = set()
    for z0 in range(p):
        c = {}
        ok = True
        for mu in (r, s, r + s, r - s):
            xp, xm = (1 + mu * z0) % p, (1 - mu * z0) % p
            a, b = chi(xp, p), chi(xm, p)
            if a and b and a != b: ok = False; break
            c[mu] = a if a else b
        if not ok: continue
        if c[r] != c[s]: continue
        if c[r + s] * c[r - s] != c[r]: continue
        img.add((c[r + s], c[r - s]))
    return img

def filter_lists(r, s, DT, DL, detail=False):
    pairs = [(a, b) for a in DT for b in DL]
    killed = {}
    for p in P3:
        if all(q in killed for q in pairs): break
        I = image_upper(r, s, p)
        for (a, b) in pairs:
            if (a, b) in killed: continue
            if (chi(a, p), chi(b, p)) not in I:
                killed[(a, b)] = (p, sorted(I))
    alive = [q for q in pairs if q not in killed]
    return alive, killed

def main():
    t0 = time.time()
    out = {"PMAX": PMAX, "SMAX": SMAX}
    # (1) три наклона и знаки
    three = {}
    for (r, s) in [(126, 451), (-126, 451), (73, 362), (-73, 362), (265, 298)]:
        DT, DL = admissible_A(r - s), admissible_A(r + s)
        alive, killed = filter_lists(r, s, DT, DL)
        I7 = sorted(image_upper(r, s, 7))
        three["%d/%d" % (r, s)] = dict(codex_pairs=[(a, b) for a in DT for b in DL], alive=alive,
                                       killed={str(k): v for k, v in killed.items()}, I_7=I7)
        print("%d/%d: пары Codex %s -> остаются %s; исключены %s; I_7=%s" % (
            r, s, [(a, b) for a in DT for b in DL], alive,
            {str(k): v[0] for k, v in killed.items()}, I7), flush=True)
    out["three"] = three
    # (2) массово
    tot = nontriv = became_trivial = still_nontriv = 0
    pairs_before = pairs_after = 0
    examples_left = []
    for s in range(2, SMAX + 1):
        for r in range(1, s):
            if gcd(r, s) != 1: continue
            tot += 1
            DT, DL = admissible_A(r - s), admissible_A(r + s)
            if DT == [1] and DL == [1]: continue
            nontriv += 1
            alive, killed = filter_lists(r, s, DT, DL)
            pairs_before += len(DT) * len(DL); pairs_after += len(alive)
            if alive == [(1, 1)]: became_trivial += 1
            else:
                still_nontriv += 1
                if len(examples_left) < 25: examples_left.append(("%d/%d" % (r, s), alive))
        if s % 50 == 0:
            print("  s=%d t=%.1fs нетрив=%d стали тривиальны=%d" % (s, time.time() - t0, nontriv, became_trivial), flush=True)
    out["mass"] = dict(slopes=tot, codex_nontrivial=nontriv, became_trivial=became_trivial,
                       still_nontrivial=still_nontriv, pairs_before=pairs_before, pairs_after=pairs_after,
                       frac_trivial_codex=(tot - nontriv) / tot, frac_trivial_refined=(tot - still_nontriv) / tot,
                       examples_left=examples_left)
    print("МАССОВО:", json.dumps(out["mass"], ensure_ascii=False), flush=True)
    print("время %.1f" % (time.time() - t0))
    json.dump(out, open("t3c_unit_filter_P%d_S%d.json" % (PMAX, SMAX), "w"), indent=1, default=str)

if __name__ == "__main__":
    main()
