# M1 (Singular genus, char 0) отдельно: аргумент h0 | rnd
import sys, time, random
R.<k,u> = QQ[]
which = sys.argv[1]
if which == 'h0':
    f = R(sage_eval(open('h0.txt').read(), locals={'k': k, 'u': u}))
else:
    random.seed(int(1))
    f = sum(ZZ(random.randint(-9, 9)) * k^i * u^j for i in range(9) for j in range(9))
t0 = time.time()
print(which, "M1 geometric genus:", Curve(f).geometric_genus(), f"({time.time()-t0:.1f}s)", flush=True)
