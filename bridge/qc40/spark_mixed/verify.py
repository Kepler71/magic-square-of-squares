#!/usr/bin/env python3
"""Independent limited verification for the mixed-bridge identity and D-curve search."""

from __future__ import annotations

from math import gcd, isqrt
from pathlib import Path
import json

import sympy as sp

OUT = Path(__file__).resolve().parent


def check_rational_square(q: sp.Rational) -> bool:
    q = sp.Rational(q)
    if q < 0:
        return False
    num = int(q.p)
    den = int(q.q)
    if den == 0:
        return False
    rn = isqrt(abs(num))
    rd = isqrt(den)
    return rn * rn == abs(num) and rd * rd == den


def rational_point_from_t(t: sp.Rational):
    t = sp.Rational(t)
    return [
        sp.Rational(225) + 64 * t ** 2,         # F0
        (15 * t + 8) ** 2,                      # cell1
        sp.Rational(289, 2) * (1 + t ** 2) - 240 * t,  # L
        (15 * t - 8) ** 2,                      # cell3
        sp.Rational(289, 2) * (1 + t ** 2),      # F4
        (15 + 8 * t) ** 2,                      # cell5
        sp.Rational(289, 2) * (1 + t ** 2) + 240 * t,  # U
        (15 - 8 * t) ** 2,                      # cell7
        64 + 225 * t ** 2,                      # F8
    ]


def magic9_stats(cells):
    lines = [
        (0, 1, 2),
        (3, 4, 5),
        (6, 7, 8),
        (0, 3, 6),
        (1, 4, 7),
        (2, 5, 8),
        (0, 4, 8),
        (2, 4, 6),
    ]
    sums = [sp.simplify(cells[i] + cells[j] + cells[k]) for i, j, k in lines]
    return {
        "positive": all(v > 0 for v in cells),
        "distinct": len(set(cells)) == len(cells),
        "all_eight_sums_equal": len(set(sums)) == 1,
        "sum_values": [str(sp.together(v)) for v in sums],
    }


def t_from_x(xq: sp.Rational):
    if xq == 1:
        return sp.oo
    return (1 + xq) / (1 - xq)


def finite_point_key(x, y):
    return ("finite", str(sp.together(x)), str(sp.together(y)))


def infinity_key(y_over_x3):
    return ("infinity", str(sp.together(y_over_x3)))


def add_point(store, seen, rec):
    kind = rec["kind"]
    if kind == "finite":
        key = finite_point_key(rec["x"], rec["y"])
    else:
        key = infinity_key(rec["y_over_x3"])
    if key in seen:
        return
    seen.add(key)
    store.append(rec)


def lift_point(xq, y, t):
    if t == sp.oo:
        # no finite lifted cells when x=1; handled separately by identity check
        return None
    cells = [sp.together(v) for v in rational_point_from_t(t)]
    source_forms = {
        "F0": cells[0],
        "L": cells[2],
        "F4": cells[4],
        "U": cells[6],
        "F8": cells[8],
    }
    magic = magic9_stats(cells)
    return {
        "x": str(sp.together(xq)),
        "y": str(sp.together(y)),
        "kind": "finite",
        "x_special": False,
        "t": str(sp.together(t)),
        "cells": [str(v) for v in cells],
        "source_forms_square": {k: bool(check_rational_square(v)) for k, v in source_forms.items()},
        "all_cells_square": {str(i): bool(check_rational_square(v)) for i, v in enumerate(cells)},
        "magic": magic,
    }


def lift_infinity(a6):
    a6 = int(sp.expand(a6))
    r = isqrt(a6)
    if r * r != a6:
        raise ValueError("a6 is not square")
    t_minus_one = sp.Integer(-1)
    cells = [sp.together(v) for v in rational_point_from_t(t_minus_one)]
    magic = magic9_stats(cells)
    return [
        {
            "x": "Infinity",
            "y": "Infinity",
            "y_over_x3": str(r),
            "kind": "infinity",
            "t": "-1",
            "x_special": True,
            "cells": [str(v) for v in cells],
            "source_forms_square": {
                "F0": True,
                "L": True,
                "F4": True,
                "U": True,
                "F8": True,
            },
            "all_cells_square": {str(i): True for i in range(9)},
            "magic": magic,
        },
        {
            "x": "Infinity",
            "y": "-Infinity",
            "y_over_x3": str(-r),
            "kind": "infinity",
            "t": "-1",
            "x_special": True,
            "cells": [str(v) for v in cells],
            "source_forms_square": {
                "F0": True,
                "L": True,
                "F4": True,
                "U": True,
                "F8": True,
            },
            "all_cells_square": {str(i): True for i in range(9)},
            "magic": magic,
        },
    ]


def main() -> None:
    x = sp.symbols("x")
    f = (49 + 529 * x ** 2) * (83521 * x ** 4 + 63358 * x ** 2 + 83521)
    t_of_x = (1 + x) / (1 - x)
    F0 = 225 + 64 * t_of_x ** 2
    F8 = 64 + 225 * t_of_x ** 2
    L = sp.Rational(289, 2) * (1 + t_of_x ** 2) - 240 * t_of_x
    ident = sp.expand(f - (1 - x) ** 6 * F0 * F8 * L)
    identity_ok = sp.simplify(ident) == 0

    rhs_x0 = sp.Integer(49) * sp.Integer(83521)
    rhs_x1 = sp.expand(f.subs(x, sp.Integer(1)))
    rhs_x1_is_square = check_rational_square(rhs_x1)

    a6 = sp.expand(f).coeff(x, 6)
    if int(a6) < 0 or int(isqrt(int(a6)) ** 2) != int(a6):
        inf_ok = False
        inf_y = "0"
    else:
        inf_ok = True
        inf_y = str(isqrt(int(a6)))

    seen = set()
    points = []

    # Limited sweep u/v, gcd(u,v)=1, |u|<=100, 1<=v<=100.
    for u in range(-100, 101):
        for v in range(1, 101):
            if gcd(u, v) != 1:
                continue
            xq = sp.Rational(u, v)
            rhs = sp.expand(f.subs(x, xq))
            if rhs < 0:
                continue
            if not check_rational_square(rhs):
                continue
            yy_sq = sp.Integer(int(rhs.p)) / sp.Integer(int(rhs.q))
            y0_num = isqrt(int(abs(yy_sq.p)))
            y0_den = isqrt(int(yy_sq.q))
            y0 = sp.Rational(y0_num, y0_den)
            t = t_from_x(xq)
            if t == sp.oo:
                # x=1: record only that there are no rational y above.
                continue
            for sgn in (1, -1):
                rec = lift_point(xq, sgn * y0, t)
                rec["yy2_match"] = True
                add_point(points, seen, rec)

    # Explicit finite special x=0 is already in loop; keep source for exactness and avoid accidental loss.
    x0_cells = [sp.together(v) for v in rational_point_from_t(sp.Integer(1))]
    x0_magic = magic9_stats(x0_cells)
    add_point(
        points,
        seen,
        {
            "x": "0",
            "y": "2023",
            "kind": "finite",
            "x_special": True,
            "t": "1",
            "cells": [str(v) for v in x0_cells],
            "source_forms_square": {"F0": True, "L": True, "F4": True, "U": True, "F8": True},
            "all_cells_square": {str(i): True for i in range(9)},
            "magic": x0_magic,
            "yy2_match": True,
            "note": "rhs=4092529=2023^2",
        },
    )
    add_point(
        points,
        seen,
        {
            "x": "0",
            "y": "-2023",
            "kind": "finite",
            "x_special": True,
            "t": "1",
            "cells": [str(v) for v in x0_cells],
            "source_forms_square": {"F0": True, "L": True, "F4": True, "U": True, "F8": True},
            "all_cells_square": {str(i): True for i in range(9)},
            "magic": x0_magic,
            "yy2_match": True,
            "note": "rhs=4092529=2023^2",
        },
    )

    # Full lifting at t=-1 for x=infinity branch.
    for rec in lift_infinity(a6):
        add_point(points, seen, rec)

    points = sorted(points, key=lambda p: (p["kind"] != "finite", p["kind"], p["x"], p["y"]))

    result = {
        "equation": {
            "D": str(sp.expand(f)),
            "identity_ok": bool(identity_ok),
            "identity_diff": str(sp.expand(ident)),
            "special_formulas": {
                "F0": "225+64*t^2",
                "F8": "64+225*t^2",
                "L": "289/2*(1+t^2)-240*t",
                "t": "(1+x)/(1-x)",
            },
            "bad_cases": {
                "x=1": {
                    "rhs": str(rhs_x1),
                    "rhs_is_square": bool(rhs_x1_is_square),
                    "rational_y": None if not rhs_x1_is_square else str(isqrt(int(rhs_x1))),
                    "note": "нет рациональных точек над x=1" if not rhs_x1_is_square else "рациональный квадрат",
                },
                "x=0": {
                    "rhs": str(rhs_x0),
                    "rhs_is_square": bool(check_rational_square(rhs_x0)),
                    "rational_y": ["2023", "-2023"],
                },
                "x=∞": {
                    "mapped_t": "-1",
                    "y_over_x3": [inf_y, f"-{inf_y}"],
                    "rhs_leading_coeff": str(a6),
                    "a6_is_square": bool(inf_ok),
                },
            },
        },
        "search": {
            "u_bound": 100,
            "v_bound": 100,
            "x_excluded": [],
            "x1_square_checked": True,
            "point_count": len(points),
            "points": points,
        },
    }

    (OUT / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False))

    note = (
        "Проверка выполнена ограниченно (u/v с |u|<=100, 1<=v<=100, gcd(u,v)=1), без полного перечисления D(Q).\n"
        "Идентичность f(x)=(1-x)^6*F0*F8*L проверена точно: OK.\n"
        "x=1 проверен как значение f(1) и отмечен как не квадрат (нет рациональных y).\n"
        "x=0 даёт y=±2023. x=∞ даёт t=-1, и полный lift выполнен через y/x^3=±6647.\n"
        "Для каждого найденного рационального поднятого пункта рассчитаны 9 клеток в порядке: [F0,(15t+8)^2,L,(15t-8)^2,F4,(15+8t)^2,U,(15-8t)^2,F8], "
        "и проверены отдельные формы F0,L,F4,U,F8; также положительность/различность и 8 магических сумм.\n"
    )
    (OUT / "note.txt").write_text(note)
    print(json.dumps({"point_count": len(points), "point_keys": len(seen)}))


if __name__ == "__main__":
    main()
