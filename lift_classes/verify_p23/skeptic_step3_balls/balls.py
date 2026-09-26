# -*- coding: utf-8 -*-
"""
Скептик, линза «Шаг 3» (p = 2, 3) теоремы Codex об остаточных классах.
Независимая реализация (не импортирует ни p23lib.py, ни файлы Codex).

Локальное утверждение, которое проверяется (не зависит от r, s):
  для b, c в Q_p с ненулевыми клетками 1 + i b + j c, если ВОСЕМЬ произведений
  по строкам / столбцам / диагоналям сетки (i,j) in {-1,0,1}^2 — квадраты в Q_p, то
     p = 2:  v2(b) >= 3 и v2(c) >= 3   (эквивалентно v2(z) >= 3 при b=rz, c=sz, (r,s)=1),
     p = 3:  v3(b) >= 1 и v3(c) >= 1   (эквивалентно v3(z) >= 1),
  и тогда все девять клеток — квадраты в Q_p.

Метод: ПОЛНЫЙ перебор шаров. Любая пара (b,c) с min(v(b),v(c)) = -m (m>0), либо m=0 и
b,c целые, записывается как b = B/p^m, c = C/p^m, B, C in Z_p (при m>0: не оба делятся на p).
Шар = (B mod p^N, C mod p^N). Для каждой клетки X = p^m + iB + jC вычисляется МНОЖЕСТВО
классов Q_p^*/Q_p^*2, которые X принимает на шаре (точно: X пробегает весь класс X0 + p^N Z_p),
затем сдвиг на класс p^{-m}. Тест совместности («ослабленный»): существует ли решение
восьми линейных уравнений классов с x_ij из этих множеств. Совместность каждого шара
проверяется по ВСЕМУ пространству решений линейной системы (оно вычислено перебором 2^9).
Ослабление (клетки независимы) может лишь ЗАВЫСИТЬ число «возможных» шаров;
поэтому «шар опровергнут» — строгий вывод для всех точек шара.

Для каждого «возможного» шара ищется ТОЧНЫЙ рациональный свидетель (b,c) в шаре
с ненулевыми клетками и восемью квадратами в Q_p (точная арифметика целых).

Все циклы — с явными конечными границами. Печать прогресса.
"""
import sys, json, time
import numpy as np

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {c: n for n, c in enumerate(CELLS)}
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)] +      # строки (i фикс.)
         [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)] +      # столбцы (j фикс.)
         [[IDX[(t, t)] for t in (-1, 0, 1)], [IDX[(t, -t)] for t in (-1, 0, 1)]])
CENTER = IDX[(0, 0)]


# ---------------- Часть A: ядро системы классов над F_2 (перебор 2^9) ----------------
def f2_kernel():
    ker = []
    for x in range(512):                       # конечная граница 2^9
        bits = [(x >> k) & 1 for k in range(9)]
        if all((bits[a] ^ bits[b] ^ bits[c]) == 0 for a, b, c in LINES):
            ker.append(bits)
    return ker


def basis_of(vecs):
    """Гауссов базис подпространства F_2^9 (вектора как списки битов)."""
    rows = []
    for v in vecs:
        w = int(''.join(map(str, v[::-1])), 2)
        for r in rows:
            w = min(w, w ^ r)
        if w:
            rows.append(w)
    return [[(r >> k) & 1 for k in range(9)] for r in rows]


# ---------------- p-адические классы ----------------
def vp_int(x, p):
    assert x != 0
    v = 0
    while x % p == 0:          # конечно: |x| убывает
        x //= p
        v += 1
    return v, x


def class_code(x, p):
    """Класс ненулевого целого x в Q_p^*/Q_p^*2, код-число.
    p=2: бит0 = v mod 2, бит1 = [u = 3 mod 4], бит2 = [u = +-3 mod 8]  (группа F_2^3)
    p=3: бит0 = v mod 2, бит1 = [u = 2 mod 3]                           (группа F_2^2)"""
    v, u = vp_int(x, p)
    if p == 2:
        return (v & 1) | ((1 if u % 4 == 3 else 0) << 1) | ((1 if u % 8 in (3, 5) else 0) << 2)
    return (v & 1) | ((1 if u % 3 == 2 else 0) << 1)


def gsize(p):
    return 8 if p == 2 else 4


def is_square_Qp_rat(num, den, p):
    """num/den != 0 — квадрат в Q_p? (точно)"""
    return class_code(num, p) ^ class_code(den, p) == 0


def residue_mask_table(p, N):
    """mask[X0] = множество (битовая маска) классов, принимаемых X на X0 + p^N Z_p."""
    M = p ** N
    tab = np.zeros(M, dtype=np.uint16)
    G = gsize(p)
    full = (1 << G) - 1
    for X0 in range(M):                         # конечная граница p^N
        if X0 == 0:
            tab[X0] = full
            continue
        w, u0 = vp_int(X0, p)                   # w < N
        known = N - w                            # u известен по модулю p^known
        if p == 2:
            if known >= 3:
                cls = {class_code(X0, 2)}
            else:
                mod = 2 ** known
                cls = set()
                for u in range(1, 8, 2):         # все нечётные u mod 8, согласованные с u0 mod 2^known
                    if (u - u0) % mod == 0:
                        cls.add(class_code((2 ** w) * u, 2))
        else:
            cls = {class_code(X0, 3)}            # known >= 1 достаточно: класс = (w mod 2, u mod 3)
        m = 0
        for c in cls:
            m |= 1 << c
        tab[X0] = m
    return tab


def shift_mask_table(tab, shift):
    """Маска классов после умножения на класс shift (XOR)."""
    if shift == 0:
        return tab.copy()
    out = np.zeros_like(tab)
    for c in range(16):
        out |= (((tab >> c) & 1) << (c ^ shift)).astype(np.uint16)
    return out


def solutions(p, kerbasis):
    """Все решения системы классов со значениями в G = F_2^g: ker (x) G."""
    G = gsize(p)
    d = len(kerbasis)
    sols = []
    for t in range(G ** d):                     # конечная граница |G|^d
        coeffs = [(t // G ** l) % G for l in range(d)]
        x = [0] * 9
        for l in range(d):
            for k in range(9):
                if kerbasis[l][k]:
                    x[k] ^= coeffs[l]
        sols.append(tuple(x))
    return sorted(set(sols))


def exact_witness(p, m, B0, C0, N, tries=40):
    """Ищем точную рациональную точку в шаре: B = B0 + p^N t, C = C0 + p^N u, |t|,|u| малы.
    Возвращаем (B, C) с ненулевыми клетками и восемью квадратами в Q_p, либо None."""
    pN = p ** N
    for t in range(-tries, tries + 1):           # конечные границы
        for u in range(-3, 4):
            B = B0 + pN * t
            C = C0 + pN * u
            X = [p ** m + i * B + j * C for (i, j) in CELLS]
            if any(x == 0 for x in X):
                continue
            ok = True
            for a, b, c in LINES:
                prod = X[a] * X[b] * X[c]            # знаменатель p^{3m}: класс p^{m}
                num = prod
                den = p ** (3 * m)
                if not is_square_Qp_rat(num, den, p):
                    ok = False
                    break
            if ok:
                return (B, C)
    return None


def run(p, m, N, sols, vtab_cache, masktab_cache, conclusion_thr):
    pN = p ** N
    key = (p, N)
    if key not in masktab_cache:
        masktab_cache[key] = residue_mask_table(p, N)
    base = masktab_cache[key]
    # класс p^{-m} = класс p^m: бит0 = m mod 2
    tab = shift_mask_table(base, m & 1)
    if key not in vtab_cache:
        vt = np.full(pN, N, dtype=np.int64)
        for x in range(1, pN):                  # конечная граница
            vt[x] = vp_int(x, p)[0]
        vtab_cache[key] = vt
    vt = vtab_cache[key]

    # центр X = p^m точно, его класс после сдвига тривиален (маска {0}) — решения с x_00 != 0
    # всё равно отвергаются маской центра; отфильтровываем их заранее ради скорости
    assert int(tab[pow(p, m, pN)] if m < N else tab[0]) >= 1
    solarr = np.array([s for s in sols if s[CENTER] == 0], dtype=np.int64)
    stats = dict(p=p, m=m, N=N, balls=0, refuted=0, possible_concl=0,
                 possible_violating=0, ambiguous=0, possible_concl_not_all9=0,
                 witness_found=0, witness_missing=0)
    violating_examples = []
    pm = p ** m
    Cs = np.arange(pN, dtype=np.int64)
    chunk = max(1, (1 << 21) // pN)
    for B_start in range(0, pN, chunk):         # конечная граница
        Bs = np.arange(B_start, min(pN, B_start + chunk), dtype=np.int64)
        BB, CC = np.meshgrid(Bs, Cs, indexing='ij')
        BB = BB.ravel(); CC = CC.ravel()
        if m > 0:
            keep = (BB % p != 0) | (CC % p != 0)
            BB = BB[keep]; CC = CC[keep]
        nb = BB.size
        if nb == 0:
            continue
        masks = np.empty((nb, 9), dtype=np.uint16)
        for k, (i, j) in enumerate(CELLS):
            X0 = (pm + i * BB + j * CC) % pN
            masks[:, k] = tab[X0]
        possible = np.zeros(nb, dtype=bool)
        for sol in solarr:                      # конечная граница
            ok = np.ones(nb, dtype=bool)
            for k in range(9):
                ok &= ((masks[:, k] >> sol[k]) & 1).astype(bool)
            possible |= ok
        vB = vt[BB % pN] - m
        vC = vt[CC % pN] - m
        # решаемость вывода на данной точности
        Bzero = (BB % pN == 0); Czero = (CC % pN == 0)
        viol = ((~Bzero) & (vB < conclusion_thr)) | ((~Czero) & (vC < conclusion_thr))
        concl = ((Bzero & (N - m >= conclusion_thr)) | ((~Bzero) & (vB >= conclusion_thr))) & \
                ((Czero & (N - m >= conclusion_thr)) | ((~Czero) & (vC >= conclusion_thr)))
        amb = ~(viol | concl)
        stats['balls'] += int(nb)
        stats['refuted'] += int((~possible).sum())
        stats['possible_concl'] += int((possible & concl).sum())
        stats['possible_violating'] += int((possible & viol).sum())
        stats['ambiguous'] += int((possible & amb).sum())
        # на шарах с выводом: все 9 клеток должны быть ОДНОЗНАЧНО квадратами (маска == {0})
        pc = possible & concl
        if pc.any():
            all9 = np.all(masks[pc] == 1, axis=1)
            stats['possible_concl_not_all9'] += int((~all9).sum())
        idx = np.nonzero(possible & viol)[0]
        for q in idx[:20]:
            violating_examples.append((int(BB[q]), int(CC[q])))
        # точные свидетели для небольшой выборки «возможных» шаров
        idxp = np.nonzero(possible)[0]
        for q in idxp[:25]:
            w = exact_witness(p, m, int(BB[q]), int(CC[q]), N)
            if w is None:
                stats['witness_missing'] += 1
            else:
                stats['witness_found'] += 1
    return stats, violating_examples


def main():
    t0 = time.time()
    out = {}
    ker = f2_kernel()
    kb = basis_of(ker)
    ker_center0 = [v for v in ker if v[CENTER] == 0]
    print(f"[A] ядро 8 линейных уравнений над F_2 (9 переменных): |ker|={len(ker)}, dim={len(kb)}")
    print(f"[A] с центром 0: {len(ker_center0)} векторов")
    for v in ker_center0:
        print("    ", [v[IDX[(i, j)]] for i in (-1, 0, 1) for j in (-1, 0, 1)])
    # сверка с формой S U S+U / U 0 U / S+U U S
    pat = []
    for S in (0, 1):
        for U in (0, 1):
            M = {(-1, -1): S, (-1, 0): U, (-1, 1): S ^ U, (0, -1): U, (0, 0): 0, (0, 1): U,
                 (1, -1): S ^ U, (1, 0): U, (1, 1): S}
            pat.append([M[c] for c in CELLS])
    okA = sorted(map(tuple, ker_center0)) == sorted(map(tuple, pat))
    print(f"[A] совпадает с матрицей Codex (S,U): {okA}")
    # T = (-1,0),(1,1),(0,-1) ; L = (-1,0),(1,-1),(0,1)
    for S in (0, 1):
        for U in (0, 1):
            M = dict(zip(CELLS, pat[2 * S + U]))
            T = M[(-1, 0)] ^ M[(1, 1)] ^ M[(0, -1)]
            L = M[(-1, 0)] ^ M[(1, -1)] ^ M[(0, 1)]
            assert T == S and L == S ^ U, (S, U, T, L)
    print("[A] [T]=S, [L]=S+U проверено на всех 4 образцах")
    out['A'] = dict(ker_size=len(ker), dim=len(kb), center0=len(ker_center0), matches_codex=okA)

    cfg = {2: dict(thr=3, Ns=range(1, 11), ms=range(0, 6)),
           3: dict(thr=1, Ns=range(1, 9), ms=range(0, 4))}
    only = [int(a) for a in sys.argv[1:]] or [2, 3]
    vcache, mcache = {}, {}
    for p in only:
        sols = solutions(p, kb)
        print(f"\n[p={p}] решений системы классов в G^9: {len(sols)}; с центром 0: "
              f"{sum(1 for s in sols if s[CENTER] == 0)}", flush=True)
        res = []
        for m in cfg[p]['ms']:
            for N in cfg[p]['Ns']:
                if p == 3 and m == 3 and N == 8:
                    continue  # бюджет времени; m=3 до N=7
                t1 = time.time()
                st, ex = run(p, m, N, sols, vcache, mcache, cfg[p]['thr'])
                st['sec'] = round(time.time() - t1, 2)
                st['violating_examples'] = ex[:5]
                res.append(st)
                print(f"  p={p} m={m} N={N}: шаров={st['balls']}, опровергнуто={st['refuted']}, "
                      f"возм.(вывод верен)={st['possible_concl']}, ВОЗМ.НАРУШАЮЩИХ={st['possible_violating']}, "
                      f"неразреш.={st['ambiguous']}, не-все-9={st['possible_concl_not_all9']}, "
                      f"свид.={st['witness_found']}/{st['witness_found'] + st['witness_missing']}, "
                      f"{st['sec']}s, всего {time.time() - t0:.0f}s", flush=True)
                if ex:
                    print("     примеры (B,C):", ex[:5], flush=True)
        out[str(p)] = res
    with open(f"balls_{'_'.join(map(str, only))}.json", "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(f"\nГотово за {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
