"""(1) Списки допустимых d_T, d_L по формулировке теоремы Codex; два независимых способа.
A: решето наименьших делителей + перебор подмножеств простых p=1 mod 4 из r-s (r+s).
B: перебор d=1 mod 24, d | m, d бесквадратное, d -- сумма двух квадратов
   (для нечётного бесквадратного d это равносильно «все простые = 1 mod 4»).
Дополнительно (НЕ утверждение Codex): фильтр Q_7 -- (d/7)=+1 (см. w3, доказательство P7).
Все циклы ограничены явно."""
import json, sys, time
from math import gcd, isqrt
from fractions import Fraction

SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
SMAX_B = 600
NMAX = 2 * SMAX + 2
t0 = time.time()

spf = list(range(NMAX + 1))
for i in range(2, isqrt(NMAX) + 1):
    if spf[i] == i:
        for j in range(i * i, NMAX + 1, i):
            if spf[j] == j:
                spf[j] = i


def primes_of(n):
    n = abs(n)
    ps = []
    while n > 1:
        p = spf[n]
        ps.append(p)
        while n % p == 0:
            n //= p
    return ps


def listA(m):
    P = [p for p in primes_of(m) if p % 4 == 1]
    out = []
    for mask in range(1 << len(P)):
        d = 1
        for k, p in enumerate(P):
            if mask >> k & 1:
                d *= p
        if d % 24 == 1:
            out.append(d)
    return sorted(out)


def sqfree(d):
    k = 2
    while k * k <= d:
        if d % (k * k) == 0:
            return False
        k += 1
    return True


def sum2sq(d):
    x = 0
    while x * x <= d:
        y2 = d - x * x
        y = isqrt(y2)
        if y * y == y2:
            return True
        x += 1
    return False


def listB(m):
    m = abs(m)
    return [d for d in range(1, m + 1, 24) if m % d == 0 and sqfree(d) and sum2sq(d)]


def q7(lst):
    return [d for d in lst if pow(d, 3, 7) == 1]


res = {"SMAX": SMAX, "SMAX_B": SMAX_B}

# --- три наклона и сверка с Codex
codex = json.load(open("/home/kep/magicKube/lift_classes/codex/proof_progress_20260926/lift/class_lists.json"))
three = {}
for key in ("126/451", "73/362", "265/298"):
    r, s = map(int, key.split("/"))
    rec = {}
    for sg in (1, -1):
        rr = sg * r
        DT, DL = listA(rr - s), listA(rr + s)
        assert DT == listB(rr - s) and DL == listB(rr + s)
        rec["%d/%d" % (rr, s)] = {"r-s": rr - s, "r+s": rr + s, "D_T": DT, "D_L": DL,
                                  "D_T_Q7": q7(DT), "D_L_Q7": q7(DL),
                                  "mod7_of_nontriv": {d: d % 7 for d in DT + DL if d > 1}}
    rec["codex_T"], rec["codex_L"] = codex[key]["T_classes"], codex[key]["L_classes"]
    rec["match_codex"] = (rec[key]["D_T"] == rec["codex_T"] and rec[key]["D_L"] == rec["codex_L"])
    three[key] = rec
    print(key, json.dumps(rec), flush=True)
res["three"] = three
res["all_match_codex"] = all(v["match_codex"] for v in three.values())

# --- массовый перебор 1 <= r < s <= SMAX, gcd=1 (r<->s не меняет списки; r->-r их меняет местами)
LEVELS = [300, 600, 1000, 2000]
st = {S: dict(slopes=0, both_triv=0, only_T_nontriv=0, only_L_nontriv=0, both_nontriv=0,
              both_triv_Q7=0) for S in LEVELS if S <= SMAX}
mismatchAB = 0; checkedB = 0; viol73 = 0; checked73 = 0; violgcd = 0; violcong = 0
maxpairs = (0, None)
ex_nontriv = []
for s in range(2, SMAX + 1):
    for r in range(1, s):
        if gcd(r, s) != 1:
            continue
        DT, DL = listA(r - s), listA(r + s)
        if s <= SMAX_B:
            checkedB += 1
            if DT != listB(r - s) or DL != listB(r + s):
                mismatchAB += 1
        for d in DT + DL:
            if not (d % 24 == 1 and d > 0 and sqfree(d)):
                violcong += 1
        for a in DT:
            for b in DL:
                if gcd(a, b) != 1:
                    violgcd += 1
        if abs(r - s) < 73:
            checked73 += 1
            viol73 += (DT != [1])
        if abs(r + s) < 73:
            checked73 += 1
            viol73 += (DL != [1])
        npairs = len(DT) * len(DL)
        if npairs > maxpairs[0]:
            maxpairs = (npairs, (r, s, DT, DL))
        if len(ex_nontriv) < 12 and (DT != [1] or DL != [1]):
            ex_nontriv.append((r, s, DT, DL))
        for S in st:
            if s <= S:
                x = st[S]
                x["slopes"] += 1
                tT, tL = DT == [1], DL == [1]
                if tT and tL:
                    x["both_triv"] += 1
                elif tL:
                    x["only_T_nontriv"] += 1
                elif tT:
                    x["only_L_nontriv"] += 1
                else:
                    x["both_nontriv"] += 1
                if q7(DT) == [1] and q7(DL) == [1]:
                    x["both_triv_Q7"] += 1
    if s % 250 == 0:
        print("  s=%d  t=%.1fs" % (s, time.time() - t0), flush=True)
for S, x in st.items():
    x["frac_both_triv"] = round(x["both_triv"] / x["slopes"], 6)
    x["frac_both_triv_Q7"] = round(x["both_triv_Q7"] / x["slopes"], 6)
    print("s<=%d:" % S, x, flush=True)
res.update(stats=st, mismatchAB=mismatchAB, checkedB=checkedB, viol73=viol73, checked73=checked73,
           violgcd=violgcd, violcong=violcong, maxpairs=maxpairs, first_nontriv=ex_nontriv)
print("A/B сверено наклонов:", checkedB, "расхождений:", mismatchAB)
print("следствие 73 (|r-s|<73 => D_T={1}; |r+s|<73 => D_L={1}): проверено", checked73, "нарушений", viol73)
print("нарушений gcd:", violgcd, " нарушений d=1 mod24/бесквадр.:", violcong)
print("макс. число пар:", maxpairs)
print("совпадение с таблицей Codex:", res["all_match_codex"])
res["time_s"] = round(time.time() - t0, 1)
print("время %.1f с" % res["time_s"])
json.dump(res, open("/home/kep/magicKube/lift_classes/recompute/v4_claude_indep/w1_lists.json", "w"), indent=1)
