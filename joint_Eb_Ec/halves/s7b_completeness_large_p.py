# -*- coding: utf-8 -*-
"""s7b: проверка полноты модели (L1)+(L2) на крупных p, где при фиксированных (b, c) много допустимых
вычетов a mod p (в s7 при p <= 211 эмпирическая оболочка для типов b+-c, b+2c, generic ограничена малым
числом допустимых a mod p, а не соотношениями). Точность p-адики снижена до 14 (достаточно: v_p(t) <= 3).
Запуск: python3 s7b_completeness_large_p.py
"""
import json, random, time
import s3_local_relations as S
from s7_joint_relation_count import counts

S.NPREC = 14
random.seed(77)
out = []
t0 = time.time()
for p, typ, npts in [(1009, "b", 600), (1009, "c", 600), (1009, "bc", 400), (2003, "b+c", 800), (2003, "b-c", 800),
                     (4001, "b+2c", 1000), (4001, "2b-c", 1000), (50021, "generic", 1200)]:
    res0 = S.find_point_type(p, typ, tries=400000)
    if res0 is None:
        print(p, typ, "точек нет", flush=True)
        continue
    b, c = res0[:2]
    emp, avals = [], set()
    for _ in range(npts):
        res = S.find_point(p, b, c, tries=20000)
        if res is None:
            continue
        a, r = res
        cl = S.classes(r, p)
        assert S.consistent(cl, S.pred_list(b, c, p, S.true_sbit(r, p)), p)
        emp.append(S.bits_of(cl))
        avals.add(a % p)
    mod = S.model_vectors(b, c, p, nsamp=4000)
    rec = {"p": p, "type": typ, "b": b, "c": c, "empirical_points": len(emp), "distinct_a_mod_p": len(avals),
           "model": counts(mod), "empirical": counts(emp), "dim_union": S.affine_rank(mod + emp)}
    out.append(rec)
    print(p, typ, "a mod p:", len(avals), "| модель:", rec["model"]["dim_total"], rec["model"]["joint_relations_all_curves"],
          "| эмпирика:", rec["empirical"]["dim_total"], rec["empirical"]["joint_relations_all_curves"],
          "| union", rec["dim_union"], "%.0fs" % (time.time() - t0), flush=True)
json.dump(out, open("s7b_completeness_large_p.json", "w"), ensure_ascii=False, indent=1)
