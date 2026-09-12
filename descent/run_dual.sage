import sys, time
load('/home/kep/magicKube/descent/run_sections.sage') if False else None
load('/home/kep/magicKube/descent/ek_descent.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
random.seed(int(11))
which = sys.argv[1:] or ['s15', 'n17', 'n61']
data = {'s15': ('s=1/5 (17,7,13)', 109, 229, 65), 'n17': ('(23,7,17)', 169, 409, 34), 'n61': ('(71,49,61)', 3061, 4381, 671)}
for key in which:
    name, A, C, D = data[key]
    t0 = time.time()
    k, rts = E1_roots(A, C, D)
    Cv = Curve3(k, rts)
    Sel = Cv.selmer_full()
    print(f"{name}: full 2-descent dim Sel2 = {Sel.dimension()} -> rank <= {Sel.dimension()-2}")
    for i in range(3):
        s1, s2 = Cv.isogeny_descent_dual(i)
        print(f"   isogeny kernel <(e{i+1},0)> (exact images + Hilbert duality): Sel(alpha) {s1}, Sel(alpha') {s2} -> rank <= {s1+s2-2}")
        if key != 'n61' or '--check' in sys.argv:
            try:
                chk = Cv.check_duality(i)
                bad = [c for c in chk if not c[4] or c[2] != c[3]]
                print(f"      duality check (sampled Im alpha' vs annihilator): {'OK' if not bad else bad}")
            except Exception as ex:
                print("      duality check failed:", ex)
    print(f"   ({time.time()-t0:.0f}s)")
