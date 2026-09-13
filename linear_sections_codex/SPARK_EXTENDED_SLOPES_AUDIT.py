#!/usr/bin/env python3
"""Certificate audit for extended_slopes_7_12.json.

Per selected `witness_index` this script validates certificate trials independently
(using exact rational arithmetic):
- chosen 3 coefficients belong to the 9-slope set,
- A = product of the 3 coefficients,
- raw ainvs define the claimed curve and torsion points lie on it,
- reconstructed p = X/A matches an all_p_tests entry and cells match
- no all-positive pairwise distinct 9-square cells in tested p values,
- eclib fields certain/upper_bound/rank,
- rank_cache[curve_key].rank_bounds == [0, 0].
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import json

SRC_JSON = Path("/home/kep/.codex/.chatgpt-projects/g-p-6aa36022d4d48191b1fc4ba89186cd9a/work/linear_sections/extended_slopes_7_12.json")
OUT_MD = Path("/home/kep/Documents/Codex/2026-09-12/spark-g1-15-8-degeneracy/outputs/extended_slopes_certificate_audit.md")


def frac(x):
    if isinstance(x, Q):
        return x
    if isinstance(x, int):
        return Q(x, 1)
    if isinstance(x, str):
        return Q(x)
    if isinstance(x, float):
        return Q(str(x))
    raise TypeError(f"Cannot parse fraction: {x} ({type(x)})")


def allowed_slopes_for_k(k: Q):
    return {
        Q(1),
        -Q(1),
        k,
        -k,
        1 + k,
        -(1 + k),
        1 - k,
        k - 1,
        Q(0),
    }


def cell_values_row_major(k: Q, p: Q):
    # row-major entries of M = [[1+p,1-p-q,1+q],[1-p+q,1,1+p-q],[1-q,1+p+q,1-p]], q = k p.
    slopes = [
        Q(1),
        -(1 + k),
        k,
        k - 1,
        Q(0),
        1 - k,
        -k,
        1 + k,
        -Q(1),
    ]
    return [1 + s * p for s in slopes]


def is_square(q: Q) -> bool:
    if q < 0:
        return False
    n, d = q.numerator, q.denominator
    return int(n ** 0.5) ** 2 == n and int(d ** 0.5) ** 2 == d


def has_all_positive_distinct_squares(vals):
    return all(v > 0 and is_square(v) for v in vals) and len(set(vals)) == len(vals)


def on_weierstrass_raw(x: Q, y: Q, ainvs):
    # y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6
    a1, a2, a3, a4, a6 = ainvs
    return y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6


def main() -> None:
    obj = json.loads(SRC_JSON.read_text())
    rank_cache = obj.get("rank_cache", {})

    # witness groups
    witness_data = defaultdict(list)
    for w in obj["slopes"]:
        wi = w["witness_index"]
        k = frac(w["k"])
        for t in w.get("tried", []):
            if "certificate" not in t:
                continue
            t = dict(t)
            t["witness_index"] = wi
            t["witness_k"] = k
            witness_data[wi].append(t)

    summary_rows = []
    all_failures = []

    for wi in sorted(witness_data):
        entries = witness_data[wi]
        passed = 0
        failed = 0

        for t in entries:
            ok = True
            msgs = []
            coeffs = [frac(c) for c in t["coefficients"]]
            k = t["witness_k"]
            allowed = allowed_slopes_for_k(k)
            for c in coeffs:
                if c not in allowed:
                    ok = False
                    msgs.append(f"coeff {c} not in 9-slope set for this k")

            A = coeffs[0] * coeffs[1] * coeffs[2]
            cert = t["certificate"]
            cert_A = frac(cert.get("A", 0))
            if cert_A != A:
                ok = False
                msgs.append(f"certificate A={cert_A} != product(coeffs)={A}")

            raw = cert.get("raw_ainvs")
            if not raw or len(raw) != 5:
                ok = False
                msgs.append("raw_ainvs missing or not length 5")
                raw_ainvs = None
            else:
                raw_ainvs = [frac(v) for v in raw]

            # curve key rank bounds
            ck = t.get("curve_key")
            rc = rank_cache.get(ck)
            if not rc:
                ok = False
                msgs.append(f"missing rank_cache entry for curve_key={ck}")
            else:
                if rc.get("rank_bounds") != [0, 0]:
                    ok = False
                    msgs.append(f"rank_bounds={rc.get('rank_bounds')} (expected [0,0])")

            # eclib at tried-level
            eclib = t.get("eclib", {})
            if eclib.get("certain") is not True:
                ok = False
                msgs.append(f"eclib.certain={eclib.get('certain')} (expected True)")
            if eclib.get("upper_bound") != 0:
                ok = False
                msgs.append(f"eclib.upper_bound={eclib.get('upper_bound')} (expected 0)")
            if eclib.get("rank") != 0:
                ok = False
                msgs.append(f"eclib.rank={eclib.get('rank')} (expected 0)")

            # torsion points and reconstructed p/cells
            all_p_tests = cert.get("all_p_tests", [])
            lookup = {}
            for row in all_p_tests:
                lookup[frac(row["p"]) ] = [frac(v) for v in row["cells"]]

            torsion_points = cert.get("torsion_points", [])
            if not raw_ainvs:
                ok = False
                msgs.append("cannot verify torsion points without raw_ainvs")
            else:
                for (xs, ys, zs) in torsion_points:
                    x = frac(xs); y = frac(ys); z = frac(zs)
                    if z == 0:
                        continue
                    if not on_weierstrass_raw(x, y, raw_ainvs):
                        ok = False
                        msgs.append(f"torsion point ({x},{y}) fails raw curve equation")
                        continue
                    if A == 0:
                        ok = False
                        msgs.append("A=0; cannot reconstruct p = X/A")
                        continue
                    p = x / A
                    tested_cells = lookup.get(p)
                    if tested_cells is None:
                        ok = False
                        msgs.append(f"no all_p_tests entry for p={p} from torsion point ({x},{y})")
                        continue

                    expected_cells = cell_values_row_major(k, p)
                    if tested_cells != expected_cells:
                        ok = False
                        msgs.append(f"cells mismatch at p={p}: expected {expected_cells}, got {tested_cells}")

            # additional check on all_p_tests list (not sign-restricted)
            if cert.get("all_p_tests"):
                for row in cert["all_p_tests"]:
                    cells = [frac(v) for v in row["cells"]]
                    if has_all_positive_distinct_squares(cells):
                        ok = False
                        msgs.append(f"all-p positive distinct squares at p={row['p']}")

            if ok:
                passed += 1
            else:
                failed += 1
                all_failures.append((wi, t, msgs))

        summary_rows.append((wi, len(entries), passed, failed))

    # write report
    lines = []
    lines.append("# Extended slopes 7_12 certificate audit")
    lines.append("")
    lines.append(f"Source: `{SRC_JSON}`")
    lines.append("")
    lines.append("## Итог по witness_index")
    for wi, total, passed, failed in summary_rows:
        lines.append(f"- `witness_index = {wi}`: проверено сертификатных попыток = `{total}`, прошло = `{passed}`, не прошло = `{failed}`")

    lines.append("")
    if all_failures:
        lines.append("## Несоответствия")
        for wi, t, msgs in all_failures:
            wkey = f"k={t['witness_k']}, witness_index={wi}, coeffs={t['coefficients']}, curve_key={t.get('curve_key')}"
            lines.append(f"- `{wkey}`")
            for m in msgs:
                lines.append(f"  - {m}")
    else:
        lines.append("## Несоответствия")
        lines.append("- нет")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
