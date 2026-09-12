import sys, itertools
sys.path.insert(0,'/home/kep/magicKube/bridge')
import refute_11_4_claude as R
from fractions import Fraction as Fr
R.PLIMIT = 200
log = R.log
# look for a curve with |Sel^2| >= 8 fully covered by explicit rational points:
# then im(E(Q)/2E(Q)) = Sel^2, so the CTP MUST vanish identically, and the test
# is no longer vacuous (non-torsion classes are involved).
found = 0
cands = []
for a in range(1, 14):
    for b in range(1, 14):
        for c in range(1, 14):
            if len({a, -b, c}) < 3: continue
            cands.append([a, -b, c])
for roots in cands:
    if found >= 2: break
    try:
        C = R.Curve(roots)
        if C.disc == 0: continue
        badp = {2}
        for t in [C.e[0]-C.e[1], C.e[0]-C.e[2], C.e[1]-C.e[2]]:
            for q in R.factorint(abs(t)): badp.add(q)
        if len(badp) > 4: continue
        S = sorted(badp)
        sel, _ = R.selmer_group(C, S, quiet=True)
        if len(sel) < 8: continue
        selset = set(tuple(x) for x in sel)
        img = set([(1,1,1)] + [tuple(t) for t in C.torsion_classes()])
        for den in range(1, 13):
            for num in range(-1500, 1501):
                x = Fr(num, den*den); v = C.f(x)
                if v <= 0: continue
                if R.sqfree(v) == 1:
                    img.add(tuple(R.sqfree(x-C.e[i]) for i in range(3)))
        img = {cl for cl in img if cl in selset}
        if len(img) != len(selset): continue
        quart = {}; worst = 1; done = 0; skip = 0
        cls = sorted(selset)
        for A in cls:
            for B in cls:
                Cc = tuple(R.sqfree(A[k]*B[k]) for k in range(3))
                gs = []; ok = True
                for cl in (A, B, Cc):
                    if cl not in quart:
                        q = R.build_quartic(C, cl, selset, verbose=False, box=800, wmax=30)
                        if q is None: ok = False; break
                        quart[cl] = q
                    gs.append(quart[cl])
                if not ok or gs[1][0] == 0: skip += 1; continue
                try:
                    v, _, _ = R.ct_pairing(C, gs[0], gs[1], gs[2], verbose=False)
                except Exception: skip += 1; continue
                done += 1
                if v == -1: worst = -1
        nont = len(selset) - 4
        log("   корни %-12s |Sel^2|=%2d (из них %d НЕторсионных) полностью покрыт точками;"
            " пар %3d, пропущено %d; есть -1: %s  %s"
            % (roots, len(selset), nont, done, skip, "ДА" if worst == -1 else "нет",
               "OK" if worst == 1 else "!!! ЛОЖНОЕ СРАБАТЫВАНИЕ"))
        found += 1
    except Exception as ex:
        continue
log("нетривиальных контрольных кривых найдено: %d" % found)
