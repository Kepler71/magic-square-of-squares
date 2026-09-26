# -*- coding: utf-8 -*-
"""s7: сколько именно СОВМЕСТНЫХ соотношений даёт каждый тип простого p (модель s3 + эмпирика).
Для аффинной оболочки W векторов классов (32 бита) и её проекций W_b, W_c, W_{b+c}, W_{b-c}
на биты линий каждой кривой: число совместных соотношений = sum dim W_curve - dim W
(соотношения, не сводящиеся к соотношениям внутри одной кривой). Отдельно -- только E_b и E_c.
Эмпирика: случайные Q_p-точки при фиксированных (p, b, c) и случайном a.
Запуск: python3 s7_joint_relation_count.py
"""
import json, random, time
from s3_local_relations import (LINES, LNAMES, find_point_type, find_point, classes, bits_of, model_vectors,
                                affine_rank, consistent, pred_list, true_sbit)

random.seed(7)
CURVE_OF = {nm: ("b" if nm.startswith("Eb_row") else "c" if nm.startswith("Ec_col")
                 else "b+c" if nm == "Ebpc_diag" else "b-c") for nm in LNAMES}


def proj(vecs, curves):
    idx = [k for k, nm in enumerate(LNAMES) for _ in range(4) if CURVE_OF[nm] in curves]
    # bits_of: для каждой линии 4 бита подряд
    idx = []
    for k, nm in enumerate(LNAMES):
        if CURVE_OF[nm] in curves:
            idx += [4 * k + s for s in range(4)]
    return [[v[i] for i in idx] for v in vecs]


def counts(vecs):
    d = affine_rank(vecs)
    per = {cv: affine_rank(proj(vecs, {cv})) for cv in ("b", "c", "b+c", "b-c")}
    dbc = affine_rank(proj(vecs, {"b", "c"}))
    return {"dim_total": d, "dim_per_curve": per,
            "joint_relations_all_curves": sum(per.values()) - d,
            "dim_Eb_Ec": dbc, "joint_relations_Eb_Ec": per["b"] + per["c"] - dbc}


if __name__ == "__main__":
  out = []
  t0 = time.time()
  for p, typ in [(3, "bc"), (5, "bc"), (13, "bc"), (13, "b"), (29, "b"), (41, "b"), (29, "c"), (41, "c"),
                 (41, "b+c"), (61, "b+c"), (97, "b+c"), (41, "b-c"), (61, "b-c"), (97, "b-c"),
                 (101, "b"), (157, "c"), (211, "b+c"), (211, "b-c"), (101, "b+2c"), (157, "b+2c"), (157, "generic")]:
      res0 = find_point_type(p, typ, tries=200000)
      if res0 is None:
          out.append({"p": p, "type": typ, "found": False})
          print(p, typ, "точек нет", flush=True)
          continue
      b, c = res0[:2]
      emp = []
      for _ in range(600):
          res = find_point(p, b, c, tries=400)
          if res is None:
              continue
          a, r = res
          cl = classes(r, p)
          assert consistent(cl, pred_list(b, c, p, true_sbit(r, p)), p)
          emp.append(bits_of(cl))
      mod = model_vectors(b, c, p, nsamp=4000)
      rec = {"p": p, "type": typ, "b": b, "c": c, "empirical_points": len(emp),
             "model": counts(mod), "empirical": counts(emp) if len(emp) > 1 else None,
             "dim_union_model_empirical": affine_rank(mod + emp)}
      out.append(rec)
      print(p, typ, "модель:", rec["model"], "| эмпирика:", rec["empirical"], "| union", rec["dim_union_model_empirical"],
            "%.0fs" % (time.time() - t0), flush=True)
  json.dump(out, open("s7_joint_relation_count.json", "w"), ensure_ascii=False, indent=1)
