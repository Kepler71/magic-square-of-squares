#!/usr/bin/env python3
"""Independent exact verification of diagonal rotation for (15,8)."""

from fractions import Fraction as Q
from pathlib import Path
import json


def F4(t: Q) -> Q:
    return Q(289, 2) * (1 + t * t)


def L(t: Q) -> Q:
    return F4(t) - Q(240) * t


def U(t: Q) -> Q:
    return F4(t) + Q(240) * t


def x_from_t(t: Q) -> Q:
    return (t - 1) / (t + 1)


def t_from_x(x: Q) -> Q:
    return (Q(1) + x) / (Q(1) - x)


def F4_x(x: Q) -> Q:
    return Q(289, 4) * (1 + x * x)


def L_x(x: Q) -> Q:
    return (Q(49) + Q(529) * x * x) / Q(4)


def U_x(x: Q) -> Q:
    return (Q(529) + Q(49) * x * x) / Q(4)


def F0_new(x: Q) -> Q:
    return Q(529) + Q(49) * x * x


def F4_new(x: Q) -> Q:
    return Q(289) * (1 + x * x)


def F8_new(x: Q) -> Q:
    return Q(49) + Q(529) * x * x


def identity_samples() -> list[dict[str, str]]:
    pts: list[Q] = [Q(1, 5), Q(2, 7), Q(-3, 4), Q(5), Q(-2, 3)]
    out = []
    for x in pts:
        if x == 1:
            continue
        t = t_from_x(x)
        if t + 1 == 0:
            continue
        d4 = F4(t) - (t + 1) ** 2 * F4_x(x)
        dL = L(t) - (t + 1) ** 2 * L_x(x)
        dU = U(t) - (t + 1) ** 2 * U_x(x)
        out.append({"x": str(x), "t": str(t), "dF4": str(d4), "dL": str(dL), "dU": str(dU)})
        assert d4 == 0 and dL == 0 and dU == 0
    return out


def c7_scaling_samples() -> list[dict[str, str]]:
    pts: list[Q] = [Q(0), Q(1, 2), Q(2, 3), Q(3), Q(-1), Q(-1, 2)]
    out = []
    for x in pts:
        out.append(
            {
                "x": str(x),
                "4*F4_x": str(Q(4) * F4_x(x)),
                "F4_new": str(F4_new(x)),
                "4*L_x": str(Q(4) * L_x(x)),
                "F8_new": str(F8_new(x)),
                "4*U_x": str(Q(4) * U_x(x)),
                "F0_new": str(F0_new(x)),
                "delta_F4": str(Q(4) * F4_x(x) - F4_new(x)),
                "delta_L": str(Q(4) * L_x(x) - F8_new(x)),
                "delta_U": str(Q(4) * U_x(x) - F0_new(x)),
            }
        )
    return out


def projective_point(t_num: Q, t_den: Q) -> dict[str, str]:
    if t_den == 0:
        return {"t": "inf", "x": "1", "comment": "x->1 when t=inf"}
    t = Q(t_num, t_den)
    if t == -1:
        return {"t": str(t), "x": "inf", "comment": "projective pole of x=(t-1)/(t+1)"}
    return {"t": str(t), "x": str(x_from_t(t)), "comment": "finite affine"}


def build_payload() -> dict:
    # projective map points
    projective = {
        "t=1": projective_point(Q(1), Q(1)),
        "t=-1": projective_point(Q(-1), Q(1)),
        "t=0": projective_point(Q(0), Q(1)),
        "t=inf": projective_point(Q(1), Q(0)),
    }

    # inverse directions (from x)
    inverse = {
        "x=-1": {"x": "-1", "t": str(t_from_x(Q(-1))), "comment": "finite regular point"},
        "x=1": {"x": "1", "t": "inf", "comment": "limit point"},
        "x=inf": {"x": "inf", "t": "-1", "comment": "inverse pole"},
    }

    old_pm = {
        "old t=+1": {
            "x": str(x_from_t(Q(1)),),
            "cells": {
                "F0'": str(F0_new(Q(0))),
                "F4'": str(F4_new(Q(0))),
                "F8'": str(F8_new(Q(0))),
            },
            "all_three_squares": True,
        },
        "old t=-1": {
            "x": "inf",
            "projective_comment": "this is pole of the x-chart (t+1=0)",
            "leading_ratio": "[49:289:529] = [7^2:17^2:23^2]",
            "all_three_squares_after_scaling": True,
        },
    }

    payload = {
        "task": "Diagonal rotation for (15,8): map x=(t-1)/(t+1)",
        "formulas_affine": {
            "F4": "289/2*(1+t^2)",
            "L": "F4-240t",
            "U": "F4+240t",
        },
        "formulas_in_x": {
            "x": "(t-1)/(t+1)",
            "t": "(1+x)/(1-x)",
            "F4/(t+1)^2": "289/4*(1+x^2)",
            "L/(t+1)^2": "(49+529x^2)/4",
            "U/(t+1)^2": "(529+49x^2)/4",
        },
        "c7_like_pair_23_7_scaled": {
            "F0'": "529+49x^2",
            "F4'": "289*(1+x^2)",
            "F8'": "49+529x^2",
            "description": "identical to (23,7)-C7 three-square part up to scale 4 and vertex relabeling",
        },
        "identity_checks": identity_samples(),
        "c7_scale_checks": c7_scaling_samples(),
        "projective_specials_t": projective,
        "projective_specials_x": inverse,
        "old_t_plus_minus_in_D": old_pm,
        "completeness_warning": "We do not claim completeness of rational points on D; only exact transformations and special-point correspondence are verified.",
    }

    return payload


if __name__ == "__main__":
    data = build_payload()
    out = Path('/home/kep/Documents/Codex/2026-09-12/spark-g1-15-8-degeneracy/outputs/rotation')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'verify_rotation.json').write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')

    md = [
        "# DIAGONAL_ROTATION_15_8",
        "",
        "Независимая проверка, только стандартный Python + Fraction.",
        "",
        "## 1) Формулы",
        "- x=(t-1)/(t+1), t=(1+x)/(1-x)",
        "- F4=289/2*(1+t^2), L=F4-240t, U=F4+240t",
        "- F4/(t+1)^2 = 289/4*(1+x^2)",
        "- L/(t+1)^2 = (49+529x^2)/4",
        "- U/(t+1)^2 = (529+49x^2)/4",
        "- Значит, после умножения на 4: 4F4'=289*(1+x^2), 4L'=49+529x^2, 4U'=529+49x^2",
        "- Это формула трёхквадратного блока для пары (23,7): F0'=529+49x^2, F4'=289*(1+x^2), F8'=49+529x^2 (с перестановкой вершин).",
        "",
        "## 2) Проверка тождеств",
        "Проверены равенства на рациональной выборке x={1/5,2/7,-3/4,5,-2/3}; все разности равны 0.",
        "",
        "## 3) Проективные особенности",
        "- t=1 -> x=0",
        "- t=0 -> x=-1",
        "- t=∞ -> x=1",
        "- t=-1 -> x=∞ (полюс карты)",
        "- x=-1 -> t=0, x=1 -> t=∞, x=∞ -> t=-1",
        "",
        "## 4) Старые t=±1 в новой карте",
        "- old t=+1 -> x=0; блок даёт F0'=529, F4'=289, F8'=49 (все квадраты).",
        "- old t=-1 -> x=∞ (проективная граница); ведущие коэффициенты 49,289,529 = 7^2,17^2,23^2 после нормировки.",
        "",
        "## 5) Важное ограничение",
        "Нельзя утверждать полноту рациональных точек D по этой проверке; здесь зафиксирована лишь корректность преобразования и соответствия особых точек.",
    ]
    (out / 'DIAGONAL_ROTATION_15_8.md').write_text("\n".join(md) + "\n", encoding='utf-8')
