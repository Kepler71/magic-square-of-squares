# Проход по чистому кручению (n = 0) для пар с квадратным центром: точки T ∈ E(Q)_tors, где f_of_point определено
# (O и y = 0 — в exceptional_fs, уже учтены). Замечание Astra: height_vectors исключал n = 0.
import sys
args = sys.argv[1:]; sys.argv = ['x']
load('/home/kep/magicKube/corners/bremner_type.sage')
for line in open(args[0]):
    b, h = map(int, line.split())
    C0, cs = conds(b, h)
    hits = []; ntors = []; err = 0
    for c1, c2, f0 in [('TL','BR',h),('TL','D',h),('BR','D',h),('BL','TR',b),('BL','D',b),('TR','D',b)]:
        try:
            bc = BaseCurve(cs[c1], cs[c2], QQ(f0))
        except Exception:
            err += 1; continue
        T = bc.E.torsion_points(); ntors.append(len(T))
        for P in T:
            f = bc.f_of_point(P)
            if f is None: continue
            assert (cs[c1][0] + cs[c1][1]*f^2).is_square() and (cs[c2][0] + cs[c2][1]*f^2).is_square()
            nsq, distinct, flags, cells = score(b, h, f, C0, cs)
            if nsq >= 7 and distinct: hits.append((c1, c2, nsq, f))
    print(f"({b},{h}) tors={ntors} err={err} hits={hits}")
