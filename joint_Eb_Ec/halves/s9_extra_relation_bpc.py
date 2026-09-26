# -*- coding: utf-8 -*-
"""s9: при p | b+c (и p | b-c) эмпирическая оболочка классов на 1 меньше модели (L1)+(L2) (s7, s7b).
Ищем лишнее соотношение: линейные функционалы (над F_2, с константой), постоянные на эмпирике, но не на модели.
Биты: для каждой линии (x, x-t)-координаты: (v mod 2, [символ = -1]).
Запуск: python3 s9_extra_relation_bpc.py [p] [type]
"""
import sys, random, json
import s3_local_relations as S

S.NPREC = 14
random.seed(99)
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2003
typ = sys.argv[2] if len(sys.argv) > 2 else "b+c"
NAMES = []
for nm in S.LNAMES:
    for k in ("x", "x-t"):
        NAMES += [nm + ":" + k + ":v", nm + ":" + k + ":chi"]


def relations(vecs):
    """Базис пространства аффинных соотношений sum c_k v_k = const, выполненных на всех векторах."""
    n = len(vecs[0])
    base = vecs[0]
    rows = []
    for v in vecs[1:]:
        rows.append([x ^ y for x, y in zip(v, base)])
    # ядро матрицы rows (над F_2): функционалы c с rows*c = 0
    M = [r[:] for r in rows]
    piv, r = [], 0
    for col in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][col]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        for i in range(len(M)):
            if i != r and M[i][col]:
                M[i] = [x ^ y for x, y in zip(M[i], M[r])]
        piv.append(col)
        r += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        c = [0] * n
        c[f] = 1
        for i, pc in enumerate(piv):
            if M[i][f]:
                c[pc] = 1
        const = sum(ci & bi for ci, bi in zip(c, base)) % 2
        basis.append((c, const))
    return basis


def holds(rel, vecs):
    c, const = rel
    return all(sum(ci & vi for ci, vi in zip(c, v)) % 2 == const for v in vecs)


res0 = S.find_point_type(p, typ, tries=400000)
b, c = res0[:2]
emp = []
for _ in range(1500):
    res = S.find_point(p, b, c, tries=20000)
    if res is None:
        continue
    a, r = res
    emp.append(S.bits_of(S.classes(r, p)))
mod = S.model_vectors(b, c, p, nsamp=6000)
R_emp = relations(emp)
extra = [rel for rel in R_emp if not holds(rel, mod)]
print("p =", p, "type", typ, "b, c =", b, c, "| символы: (-1/p) =", S.legendre(-1, p), "(2/p) =", S.legendre(2, p),
      "(b/p) =", S.legendre(b, p))
print("соотношений на эмпирике:", len(R_emp), "; из них не следуют из модели:", len(extra))
# упростим: приведём лишние соотношения по модулю модельных
R_mod = relations(mod)
for c_, const in extra:
    supp = [NAMES[k] for k, x in enumerate(c_) if x]
    print("лишнее:", " + ".join(supp), "=", const)
json.dump({"p": p, "type": typ, "b": b, "c": c, "n_emp_rel": len(R_emp), "n_extra": len(extra),
           "extra": [[NAMES[k] for k, x in enumerate(c_) if x] + ["=", const] for c_, const in extra]},
          open("s9_extra_relation_%s_%d.json" % (typ.replace("+", "p").replace("-", "m"), p), "w"),
          ensure_ascii=False, indent=1)
